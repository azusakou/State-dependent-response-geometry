#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
pa<-fread(file.path(out,'response_predictions_by_axis_primary.tsv.gz'))
u<-fread(file.path(out,'analysis_units_partition1.tsv.gz'))
pa[,squared_error:=(prediction-observed)^2]
pa[,axis_group:=ifelse(axis%in%c('apoptosis','proliferation','stress'),'core3','hallmark9')]
grouped<-pa[,.(mse=mean(squared_error)),by=.(dataset,unit_id,context,perturbation,model,axis_group)]
all12<-pa[,.(mse=mean(squared_error),axis_group='all12'),by=.(dataset,unit_id,context,perturbation,model)]
grouped<-rbindlist(list(grouped,all12),use.names=TRUE)
gm<-grouped[,.(n_units=.N,rmse=sqrt(mean(mse))),by=.(dataset,axis_group,model)]
fwrite(gm,file.path(out,'primary_response_rmse_by_axis_group.tsv'),sep='\t')
unit_meta<-u[,.(unit_id,biological_context)]
grouped<-merge(grouped,unit_meta,by='unit_id')
rp<-fread(file.path(out,'leave_study_out_predictions.tsv.gz'))[partition_id==1&task=='reliability']
set.seed(202611063L);results<-list()
comparison<-function(z,study,task,axis_group,model,reference,err_model,err_reference){
 definitions<-if(study=='Norman')c('complete_condition','deterministic_anchor_gene')else'complete_condition'
 for(definition in definitions){
 key<-z$perturbation
 if(definition=='deterministic_anchor_gene')key<-vapply(strsplit(key,'_',fixed=TRUE),function(g)sort(g[!grepl('^NegCtrl',g)])[1],character(1))
 ii<-match(key,unique(key));context<-if(study=='MIX-Seq')z$biological_context else z$context
 cc<-match(context,unique(context));ni<-max(ii);nc<-max(cc)
 for(method in c('intervention','context','crossed_intervention_context')){
  stat<-function(w)sqrt(sum(err_model*w)/sum(w))-sqrt(sum(err_reference*w)/sum(w))
  bo<-replicate(2000,{
   wi<-if(method%in%c('intervention','crossed_intervention_context'))tabulate(sample.int(ni,ni,replace=TRUE),ni)[ii]else rep(1,length(ii))
   wc<-if(method%in%c('context','crossed_intervention_context'))tabulate(sample.int(nc,nc,replace=TRUE),nc)[cc]else rep(1,length(cc))
   w<-wi*wc;if(!any(w>0))NA_real_ else stat(w)
  });bo<-bo[is.finite(bo)]
  results[[length(results)+1L]]<<-data.table(dataset=study,task=task,axis_group=axis_group,model=model,reference=reference,rmse_difference=stat(rep(1,nrow(z))),ci_low=quantile(bo,.025),ci_high=quantile(bo,.975),bootstrap=method,n_units=nrow(z),intervention_clusters=ni,intervention_definition=definition,context_clusters=nc,context_definition=if(study=='MIX-Seq')'cell_line'else'technical_GEM',valid_draws=length(bo),scope='conditional fixed OOF predictions and split; no representation/model refit; technical context is not independent biological replication')
 }
 }
}
for(s in unique(grouped$dataset))for(ag in c('core3','hallmark9','all12')){
 z<-merge(grouped[dataset==s&axis_group==ag&model=='pilot_plus_state'],grouped[dataset==s&axis_group==ag&model=='pilot_response',.(unit_id,ref_mse=mse)],by='unit_id')
 comparison(z,s,'response',ag,'pilot_plus_state','pilot_response',z$mse,z$ref_mse)
}
for(s in unique(rp$dataset)){
 z<-merge(rp[dataset==s&model=='sampling_plus_state'],rp[dataset==s&model=='sampling_aware',.(unit_id,ref_error=log_squared_error)],by='unit_id')
 comparison(z,s,'reliability_log_disagreement','all12','sampling_plus_state','sampling_aware',z$log_squared_error,z$ref_error)
}
r<-rbindlist(results);fwrite(r,file.path(out,'direct_state_value_cluster_comparisons.tsv'),sep='\t')
print(gm[axis_group=='core3'&model%in%c('zero','training_mean','pilot_identity','control_state','pilot_response','pilot_plus_state')])
print(r[bootstrap=='crossed_intervention_context'])
