#!/usr/bin/env Rscript
# Conditional analyses of archived out-of-fold predictions. No refitting is implied.
suppressPackageStartupMessages(library(data.table))
set.seed(191103)
args <- commandArgs(trailingOnly=TRUE)
root <- '.'
out <- 'results/prediction_checks'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd <- function(p) fread(file.path(root,p))
wr <- function(x,p) fwrite(x,file.path(out,p),sep='\t')
p <- rd('results/mixseq_control_models/predictions.tsv')
d <- as.data.table(readRDS('inputs/mixseq_modeling.rds'))
reg <- rd('results/task_registry/state_transition_tasks.tsv')
d <- d[pair_id %in% reg[dataset_id=='mixseq_2020'&status=='eligible',task_id]]
zcols <- grep('^z[0-9]+_mean$',names(d),value=TRUE)
# Verify the row ordering used in the archived split manifest.
base <- unique(p[,.(row_id,pair_id,cell_line)])
stopifnot(all(d$pair_id[base$row_id]==base$pair_id))
d[,row_id:=.I]
wr(d[,c('row_id','pair_id','cell_line','perturbation','time_hours',zcols,'p_apoptosis','p_proliferation','p_stress'),with=FALSE],'mixseq_benchmark_input.tsv')
# Hard applicability mask reconstructed using training-fold predictors only.
splits <- rd('results/mixseq_control_models/split_assignments.tsv')
sup <- list()
for(sp in unique(splits$split)) for(fo in unique(splits[split==sp,fold])) {
 te <- splits[split==sp&fold==fo,row_id]; tr <- setdiff(d$row_id,te)
 z <- as.matrix(d[,..zcols]); mu <- colMeans(z[tr,,drop=FALSE]); ss <- apply(z[tr,,drop=FALSE],2,sd);ss[!is.finite(ss)|ss<1e-8]<-1
 z <- sweep(sweep(z,2,mu,'-'),2,ss,'/')
 dd <- as.matrix(dist(z[tr,,drop=FALSE]));diag(dd)<-Inf
 cut <- quantile(apply(dd,1,min),.95)
 near <- vapply(te,function(i)min(sqrt(rowSums((z[tr,,drop=FALSE]-rep(z[i,],each=length(tr)))^2))),numeric(1))
 npert <- vapply(te,function(i)sum(d$perturbation[tr]==d$perturbation[i]),integer(1))
 sup[[length(sup)+1L]] <- data.table(split=sp,fold=fo,row_id=te,unseen_intervention=npert==0,source_outside_support=near>cut,kernel_support_n=npert,source_distance=near,source_cutoff=as.numeric(cut))
}
p <- merge(p,rbindlist(sup),by=c('split','fold','row_id'))
p[,hard_support:=!unseen_intervention&!source_outside_support&(method!='kernel_transition'|kernel_support_n>=10)]
resp<-c('local_response_ridge','kernel_transition','conditional_rate_network')
p[method%in%resp,reconstructed_accept:=hard_support&is.finite(uncertainty)&uncertainty<=.15]
stopifnot(all(p[method%in%resp,reconstructed_accept==!abstain]))
wr(unique(p[method%in%resp,.(split,method,row_id,pair_id,cell_line,perturbation,unseen_intervention,source_outside_support,kernel_support_n,hard_support,uncertainty,accepted=!abstain)]),'mixseq_applicability.tsv')
ends<-c('apoptosis','proliferation','stress')
for(e in ends) p[,(paste0('se_',e)):=(get(paste0('predicted_',e))-get(paste0('observed_',e)))^2]
p[,se_overall:=rowMeans(.SD),.SDcols=paste0('se_',ends)]
# Cell-line cluster bootstrap keeps all exposure/drug conditions within each sampled line.
B<-2000L
bootmetric<-function(dt,cols,fn=function(x)sqrt(mean(x)),cluster='cell_line'){
 gs<-split(seq_len(nrow(dt)),dt[[cluster]]);ng<-length(gs)
 if(!ng)return(rep(NA_real_,3))
 point<-fn(as.matrix(dt[,..cols])); if(ng<2)return(c(point,NA_real_,NA_real_))
 v<-replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);fn(as.matrix(dt[ii,..cols]))})
 c(point,unname(quantile(v,c(.025,.975),na.rm=TRUE)))
}
metrics<-list()
for(sp in unique(p$split)) for(me in unique(p$method)) {
 a<-p[split==sp&method==me]; b<-a[abstain==FALSE]
 for(e in c('overall',ends)) {
  m<-if(nrow(b))bootmetric(b,paste0('se_',e)) else rep(NA_real_,3)
  metrics[[length(metrics)+1L]]<-data.table(split=sp,method=me,endpoint=e,n_total=nrow(a),n_accepted=nrow(b),n_cell_lines=uniqueN(b$cell_line),coverage=nrow(b)/nrow(a),rmse=m[1],ci_low=m[2],ci_high=m[3],bootstrap_unit='cell_line',bootstrap_replicates=B)
 }
}
wr(rbindlist(metrics),'mixseq_metrics_ci.tsv')
# Pairwise identical accepted cases: both methods must accept a case.
paired<-list(); meth<-unique(p$method)
for(sp in unique(p$split)) for(i in 1:(length(meth)-1))for(j in (i+1):length(meth)){
 a<-p[split==sp&method==meth[i]&!abstain];b<-p[split==sp&method==meth[j]&!abstain]
 a<-merge(a,b[,.(row_id,se_other=se_overall)],by='row_id');if(!nrow(a))next
 point<-sqrt(mean(a$se_overall))-sqrt(mean(a$se_other)); gs<-split(seq_len(nrow(a)),a$cell_line);ng<-length(gs)
 vv<-if(ng>1)replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);sqrt(mean(a$se_overall[ii]))-sqrt(mean(a$se_other[ii]))})else NA_real_
 ci<-if(ng>1)quantile(vv,c(.025,.975))else c(NA_real_,NA_real_)
 paired[[length(paired)+1L]]<-data.table(split=sp,method_a=meth[i],method_b=meth[j],n_matched=nrow(a),n_cell_lines=ng,rmse_a=sqrt(mean(a$se_overall)),rmse_b=sqrt(mean(a$se_other)),difference_a_minus_b=point,ci_low=ci[1],ci_high=ci[2])
}
wr(rbindlist(paired),'mixseq_matched_comparisons.tsv')
# All six methods on their common accepted intersection (may be empty).
common<-list()
for(sp in unique(p$split)){
 ids<-p[split==sp&!abstain,.(n_methods=uniqueN(method)),by=row_id][n_methods==length(meth),row_id]
 for(me in meth){a<-p[split==sp&method==me&row_id%in%ids];m<-if(nrow(a))bootmetric(a,'se_overall')else rep(NA_real_,3)
 common[[length(common)+1L]]<-data.table(split=sp,method=me,n=nrow(a),coverage=nrow(a)/uniqueN(p[split==sp,row_id]),rmse=m[1],ci_low=m[2],ci_high=m[3])}
}
wr(rbindlist(common),'mixseq_all_model_common_cases.tsv')
curves<-list();cal<-list()
for(sp in unique(p$split))for(me in resp){
 a<-p[split==sp&method==me]; valid<-a[hard_support&is.finite(uncertainty)]
 taus<-sort(unique(c(0,.05,.10,.15,.20,.25,.30,.50,Inf,valid$uncertainty)))
 for(tau in taus){b<-valid[uncertainty<=tau];ids<-b$row_id
 for(cmp in c(me,'elastic_net','ridge_static','training_mean')){
 v<-p[split==sp&method==cmp&row_id%in%ids]
 for(e in c('overall',ends))curves[[length(curves)+1L]]<-data.table(split=sp,selection_model=me,prediction_model=cmp,threshold=tau,endpoint=e,n=nrow(v),coverage=nrow(v)/nrow(a),rmse=if(nrow(v))sqrt(mean(v[[paste0('se_',e)]]))else NA_real_)
 }}
 cal[[length(cal)+1L]]<-data.table(split=sp,method=me,n_hard_supported=nrow(valid),uncertainty_error_spearman=if(nrow(valid)>2&&sd(valid$uncertainty)>0)cor(valid$uncertainty,sqrt(valid$se_overall),method='spearman')else NA_real_)
}
wr(rbindlist(curves),'mixseq_risk_coverage.tsv');wr(rbindlist(cal),'mixseq_uncertainty_error.tsv')
# Retained versus excluded task profiles, without new eligibility claims.
allreg<-reg[dataset_id%in%c('mixseq_2020','norman_2019_combinatorial','replogle_2022_k562_essential')]
wr(allreg[,.(n=.N,median_target_cells=as.numeric(median(target_cell_count,na.rm=TRUE)),median_response_distance=median(target_state_distance,na.rm=TRUE),q25_distance=quantile(target_state_distance,.25,na.rm=TRUE),q75_distance=quantile(target_state_distance,.75,na.rm=TRUE),median_uncertainty=median(uncertainty,na.rm=TRUE)),by=.(dataset_id,status)],'retained_excluded_characteristics.tsv')
# Norman conditional OOF evaluation. These intervals do not estimate independent biological replication.
n<-rd('results/grouped_validation/norman_predictions.tsv');nb<-rd('results/extended_baselines/norman_baseline_predictions.tsv')
np<-rbindlist(list(n[,.(split,fold,held_out,combination,method,predicted,observed,observed_class,observed_binary)],nb[,.(split,fold,held_out,combination,method,predicted,observed,observed_class,observed_binary)]))
auc<-function(y,s){if(!any(y)||all(y))return(NA_real_);(sum(rank(s)[y])-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(!y))}
ap<-function(y,s){if(!any(y))return(NA_real_);o<-order(-s); yy<-y[o]; ss<-s[o];ends<-which(!duplicated(ss,fromLast=TRUE)); sum(cumsum(yy)[ends]/ends * c(diff(c(0,cumsum(yy)[ends]))))/sum(yy)}
stat<-function(a)c(rmse=sqrt(mean((a$predicted-a$observed)^2)),spearman=if(sd(a$predicted)>0&&sd(a$observed)>0)cor(a$predicted,a$observed,method='spearman')else NA_real_,auroc=auc(a$observed_binary,a$predicted),average_precision=ap(a$observed_binary,a$predicted))
ns<-list()
for(sp in unique(np$split))for(me in unique(np$method)){
 a<-np[split==sp&method==me]; point<-stat(a)
 gs<-split(seq_len(nrow(a)),a$held_out);ng<-length(gs)
 bo<-replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);stat(a[ii])})
 for(k in seq_along(point)){
 v<-bo[k,];v<-v[is.finite(v)];ci<-if(length(v))quantile(v,c(.025,.975))else c(NA_real_,NA_real_)
 ns[[length(ns)+1L]]<-data.table(split=sp,method=me,n=nrow(a),positives=sum(a$observed_binary),prevalence=mean(a$observed_binary),metric=names(point)[k],estimate=point[k],ci_low=ci[1],ci_high=ci[2],bootstrap_unit=if(sp=='leave_gene')'assigned_anchor_gene'else'primary_Pfam_family',n_clusters=ng,valid_bootstraps=length(v),bootstrap_replicates=B)
 }
}
wr(rbindlist(ns),'norman_prediction_ci.tsv')
wr(unique(np[,.(split,fold,held_out,combination,observed_binary)])[,.(n_pairs=.N,positive_pairs=sum(observed_binary),prevalence=mean(observed_binary)),by=.(split,fold)],'norman_fold_class_counts.tsv')
# Conditional class / magnitude / precision associations.
co<-rd('results/norman_control/combination_interactions.tsv')
co[,single_response_sum:=single_a_movement+single_b_movement]
co[,magnitude_quartile:=as.integer(cut(rank(single_response_sum,ties.method='first'),breaks=quantile(rank(single_response_sum,ties.method='first'),seq(0,1,.25)),include.lowest=TRUE))]
co[,precision_quartile:=as.integer(cut(rank(replicate_projection_se,ties.method='first'),breaks=quantile(rank(replicate_projection_se,ties.method='first'),seq(0,1,.25)),include.lowest=TRUE))]
wr(co[,.(n=.N,positives=sum(interaction_class=='synergistic'),antagonistic=sum(interaction_class=='antagonistic'),positive_fraction=mean(interaction_class=='synergistic'),median_single_magnitude=median(single_response_sum),median_projection_se=median(replicate_projection_se)),by=magnitude_quartile],'norman_class_by_magnitude.tsv')
wr(co[,.(n=.N,positives=sum(interaction_class=='synergistic'),antagonistic=sum(interaction_class=='antagonistic'),positive_fraction=mean(interaction_class=='synergistic'),median_single_magnitude=median(single_response_sum),median_projection_se=median(replicate_projection_se)),by=precision_quartile],'norman_class_by_precision.tsv')
wr(data.table(comparison=c('positive_class_vs_single_magnitude','positive_class_vs_projection_SE','magnitude_vs_projection_SE'),spearman=c(cor(as.numeric(co$interaction_class=='synergistic'),co$single_response_sum,method='spearman'),cor(as.numeric(co$interaction_class=='synergistic'),co$replicate_projection_se,method='spearman'),cor(co$single_response_sum,co$replicate_projection_se,method='spearman'))),'norman_class_correlations.tsv')
# Same-drug continued-exposure endpoints, retrospective support labels.
sq<-rd('results/control_atlas/observed_support_sequences.tsv')
sq[,max6:=pmax(p_apoptosis_6,p_proliferation_6,p_stress_6)]
sq[,max24:=pmax(p_apoptosis_24,p_proliferation_24,p_stress_24)]
sq[,early_target:=ends[max.col(as.matrix(.SD),ties.method='first')],.SDcols=paste0('p_',ends,'_6')]
sq[,late_target:=ends[max.col(as.matrix(.SD),ties.method='first')],.SDcols=paste0('p_',ends,'_24')]
sq[,`:=`(target_changed=early_target!=late_target,reported_positive=observed_support_policy=='continue_same_intervention_to_24h',both_precise=uncertainty_6<=.10&uncertainty_24<=.10)]
sq[,fixed_early_target_change:=vapply(seq_len(.N),function(i)sq[[paste0('p_',sq$early_target[i],'_24')]][i]-sq[[paste0('p_',sq$early_target[i],'_6')]][i],numeric(1))]
wr(sq,'continued_exposure_reassessment.tsv')
wr(sq[,.(n=.N,reported_positive=sum(reported_positive),retrospective_nonpositive=sum(!reported_positive),both_precise=sum(both_precise),max_probability_increase=sum(max24>max6),early_target_increase=sum(fixed_early_target_change>0),target_changed=sum(target_changed),median_max_change=median(max24-max6),median_fixed_target_change=median(fixed_early_target_change)),by=perturbation],'continued_exposure_by_drug.tsv')

input_dir<-out;output_pdf<-file.path(out,"selective_risk_coverage.pdf")
c <- fread(file.path(input_dir,'mixseq_risk_coverage.tsv'))
m <- fread(file.path(input_dir,'mixseq_metrics_ci.tsv'))
models <- c('local_response_ridge','kernel_transition','conditional_rate_network')
labels <- c('Local ridge','Kernel transition','Conditional rate')
colors <- c('#0072B2','#D55E00','#009E73')
axes <- c('leave_cell_line','leave_perturbation','leave_source_state','leave_target_state')
endpoints <- c('overall','apoptosis','proliferation','stress')
pdf(output_pdf,width=11,height=8.5,onefile=TRUE,family='Helvetica')
for(ax in axes){
 par(mfrow=c(2,2),mar=c(4.0,4.1,2.4,1.0),oma=c(0,0,2.0,0))
 for(ep in endpoints){
  d<-c[split==ax&endpoint==ep&n>0&prediction_model%in%c(models,'elastic_net')]
  all_elastic<-m[split==ax&endpoint==ep&method=='elastic_net']
  ymax<-max(c(.35,d$rmse,all_elastic$rmse),na.rm=TRUE)*1.08
  plot(NA,xlim=c(0,1),ylim=c(0,ymax),xlab='Coverage (all held-out conditions)',ylab='RMSE',main=if(ep=='overall')'All three endpoints'else paste0(toupper(substr(ep,1,1)),substr(ep,2,nchar(ep))))
  for(j in seq_along(models)){
   sel<-models[j]
   for(pred in c(sel,'elastic_net')){
    v<-d[selection_model==sel&prediction_model==pred][order(coverage,threshold)]
    lines(v$coverage,v$rmse,col=colors[j],lwd=1.7,lty=if(pred==sel)1 else 2)
    dot<-v[abs(threshold-.15)<1e-10]
    points(dot$coverage,dot$rmse,col=colors[j],pch=if(pred==sel)19 else 1,cex=.85)
   }
  }
  points(all_elastic$coverage,all_elastic$rmse,pch=17,col='black',cex=.9)
  if(!nrow(d))text(.46,ymax*.45,'No supported response-model predictions',cex=.8)
  legend('topright',legend=c(labels,'Dashed: matched elastic net','Elastic net, all cases'),col=c(colors,'#555555','black'),lty=c(1,1,1,2,NA),pch=c(NA,NA,NA,NA,17),lwd=c(1.7,1.7,1.7,1.7,NA),cex=.68,bty='n',y.intersp=1.05)
 }
 mtext(gsub('_',' ',ax),outer=TRUE,font=2,cex=1.15)
}
dev.off()
cat('Four pages: cell line, intervention, source state, target state.\n')

