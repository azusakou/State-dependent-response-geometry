#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
for(e in parse(file.path(out,'fit_transfer.R'))){
 if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e)
}
lambdas<-c(.01,.1,1,10,100)
d<-fread(file.path(out,'analysis_units_partition1.tsv.gz'))
axes<-unique(fread(file.path(out,'final_axis_membership.tsv'))$program)
d[,`:=`(log_n_treated=log(n_treated_A+n_treated_B),log_n_control=log(n_control_A+n_control_B))]
for(a in axes)d[,(paste0('log_C_var_',a)):=log(get(paste0('C_var_',a))+1e-4)]
features<-c(paste0('response_',axes,'_A'),paste0('C_mean_',axes),paste0('log_C_var_',axes),'log_n_treated','log_n_control','C_mean_log_library')
target<-paste0('response_',axes,'_B')
checks<-list()
for(s in unique(d$dataset)){
 a<-d[dataset!=s];b<-d[dataset==s];changed<-copy(b)
 changed[,(target):=lapply(.SD,function(x)-100*x+25),.SDcols=target]
 original<-select_predict(a,b,features,target);mutated<-select_predict(a,changed,features,target)
 delta<-max(abs(original$prediction-mutated$prediction))
 stopifnot(delta==0,identical(original$lambda,mutated$lambda),!any(a$unit_id%in%b$unit_id))
 checks[[length(checks)+1L]]<-data.table(held_out_study=s,held_out_B_mutation_max_prediction_difference=delta,selected_lambda_unchanged=TRUE,training_test_unit_disjoint=TRUE)
}
membership<-fread(file.path(out,'primary_cell_partition_membership.tsv.gz'))
stopifnot(!anyDuplicated(membership[,.(dataset,source_file,cell_id)]))
fwrite(rbindlist(checks),file.path(out,'no_label_leakage_checks.tsv'),sep='\t')
cat('Held-out B mutation has zero effect on predictions/tuning; primary cell partitions are disjoint.\n')
