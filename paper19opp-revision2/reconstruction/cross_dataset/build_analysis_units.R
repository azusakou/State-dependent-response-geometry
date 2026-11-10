#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
d<-readRDS(file.path(out,'common_cell_scores.rds'));d<-as.data.table(d)
# Chemical MIX-Seq is a distinct experimental source; its ancillary GPX4 screen
# is not relabelled as an additional independent source.
d<-d[dataset!='MIX-Seq'|task=='chemical']
sc<-grep('^score_',names(d),value=TRUE);axes<-sub('^score_','',sc)
setorder(d,dataset,context,perturbation,cell_id)
controls<-d[is_control==TRUE];treated<-d[is_control==FALSE]
counts<-treated[,.(treated_n=.N),by=.(dataset,context,biological_context,perturbation,mode)]
cn<-controls[,.(control_n=.N),by=.(dataset,context)]
counts<-merge(counts,cn,by=c('dataset','context'),all.x=TRUE)
fwrite(counts,file.path(out,'unit_eligibility.tsv'),sep='\t')
ok<-counts[treated_n>=10&control_n>=60]
stopifnot(all(is.finite(ok$control_n)))
fwrite(ok[,.(eligible_units=.N,interventions=uniqueN(perturbation),contexts=uniqueN(context),biological_contexts=uniqueN(biological_context)),by=dataset],file.path(out,'analysis_denominators.tsv'),sep='\t')
protocol<-c('Study count: three (Replogle CRISPRi K562, Norman CRISPRa K562, MIX-Seq chemical).',
'Feature universe and program genes fixed from five-study measured gene intersection; no outcome-dependent feature selection.',
'Twelve fixed RNA axes: three short signed modules and nine Hallmark positive-mean log-expression scores; overlapping axes are not independent replicates.',
'Fixed original source cell QC; no response/magnitude/agreement eligibility gate.',
'Primary eligibility: at least ten treated cells per intervention/context and sixty controls.',
'Contexts are technical GEM partitions in Replogle/Norman and experiment/horizon/cell line in MIX-Seq.',
'Controls randomly split into three nonoverlapping sets: A/B anchor responses; C independently supplies context-state predictors.',
'Treated cells randomly split into two nonoverlapping sets A/B. No post-treatment feature is called baseline.',
'Ten fixed-seed cell partitions: first is primary, others assess partition sensitivity.',
'Strict response models use C context scores alone or C context plus explicitly labelled sampling-support covariates.',
'Pilot response model A-to-B is separate from unmeasured-intervention prediction.',
'Reliability endpoint is mean squared A/B response disagreement; candidate features distinguish counts, C state/variance and pilot response.',
'Three cross-program coupling tests per study: stress(A)-proliferation(B), stress(A)-apoptosis(B), apoptosis(A)-proliferation(B).',
'Coupling adjustment uses context fixed effects, log cell counts, log common-library totals; intervention-cluster uncertainty; BH across nine primary comparisons.',
'C-context scores are contemporaneous matched-control state, not individual pretreatment state; study transfer does not establish a causal baseline effect.')
writeLines(protocol,file.path(out,'locked_protocol.txt'))
partition_count<-if(length(args)>=3L)as.integer(args[3])else 10L
for(rep in seq_len(partition_count)){
 set.seed(202611060L+rep)
 cc<-controls[paste(dataset,context)%in%paste(ok$dataset,ok$context)]
 cc[,partition:=sample(rep(LETTERS[1:3],length.out=.N)),by=.(dataset,context)]
 cs<-cc[,{
  mm<-lapply(.SD,mean);vv<-lapply(.SD,var)
  c(list(n=.N,mean_log_library=mean(log1p(common_library))),setNames(mm,paste0('mean_',axes)),setNames(vv,paste0('var_',axes)))
 },by=.(dataset,context,partition),.SDcols=sc]
 tt<-merge(treated,ok[,.(dataset,context,perturbation)],by=c('dataset','context','perturbation'))
 tt[,partition:=sample(rep(c('A','B'),length.out=.N)),by=.(dataset,context,perturbation)]
 if(rep==1L){
  membership<-rbindlist(list(cc[,.(dataset,context,perturbation,cell_id,source_file,partition,role='control')],tt[,.(dataset,context,perturbation,cell_id,source_file,partition,role='treated')]))
  stopifnot(!anyDuplicated(membership[,.(dataset,source_file,cell_id)]))
  fwrite(membership,file.path(out,'primary_cell_partition_membership.tsv.gz'),sep='\t')
 }
 ts<-tt[,c(list(n=.N,mean_log_library=mean(log1p(common_library))),setNames(lapply(.SD,mean),paste0('mean_',axes)),setNames(lapply(.SD,var),paste0('var_',axes))),by=.(dataset,context,biological_context,perturbation,mode,partition),.SDcols=sc]
 rr<-merge(ts,cs[partition%in%c('A','B')],by=c('dataset','context','partition'),suffixes=c('_treated','_control'))
 for(a in axes)rr[,(paste0('response_',a)):=get(paste0('mean_',a,'_treated'))-get(paste0('mean_',a,'_control'))]
 keepcols<-c('dataset','context','biological_context','perturbation','mode','partition','n_treated','n_control','mean_log_library_treated','mean_log_library_control',paste0('response_',axes),paste0('var_',axes,'_treated'),paste0('var_',axes,'_control'))
 rwide<-dcast(rr[,..keepcols],dataset+context+biological_context+perturbation+mode~partition,value.var=setdiff(keepcols,c('dataset','context','biological_context','perturbation','mode','partition')))
 cfeat<-cs[partition=='C'];cfeat[,partition:=NULL];setnames(cfeat,setdiff(names(cfeat),c('dataset','context')),paste0('C_',setdiff(names(cfeat),c('dataset','context'))))
 z<-merge(rwide,cfeat,by=c('dataset','context'))
 z[,`:=`(partition_id=rep,unit_id=paste(dataset,context,perturbation,sep='|'))]
 for(a in axes){
  z[,(paste0('squared_disagreement_',a)):=(get(paste0('response_',a,'_A'))-get(paste0('response_',a,'_B')))^2]
  z[,(paste0('analytic_sampling_',a)):=(get(paste0('var_',a,'_treated_A'))/n_treated_A+get(paste0('var_',a,'_treated_A'))/n_treated_B+get(paste0('C_var_',a))/n_control_A+get(paste0('C_var_',a))/n_control_B)]
 }
 z[,disagreement:=rowMeans(as.matrix(.SD)),.SDcols=paste0('squared_disagreement_',axes)]
 z[,analytic_sampling:=rowMeans(as.matrix(.SD)),.SDcols=paste0('analytic_sampling_',axes)]
 fwrite(z,file.path(out,paste0('analysis_units_partition',rep,'.tsv.gz')),sep='\t')
 if(rep==1)fwrite(cs,file.path(out,'independent_control_parts.tsv'),sep='\t')
 cat('Partition',rep,'units',nrow(z),'\n')
}
