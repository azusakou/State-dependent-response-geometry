#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
axis_spec<-fread(file.path(out,'final_axis_membership.tsv'));axes<-unique(axis_spec$program);lambdas<-c(.01,.1,1,10,100)
ridge<-function(x,y,xx,w,lambda){
 x<-as.matrix(x);xx<-as.matrix(xx);y<-as.matrix(y);w<-w/sum(w)
 mu<-colSums(x*w);ss<-sqrt(colSums(sweep(x,2,mu,'-')^2*w));ss[!is.finite(ss)|ss<1e-8]<-1
 x<-sweep(sweep(x,2,mu,'-'),2,ss,'/');xx<-sweep(sweep(xx,2,mu,'-'),2,ss,'/')
 ym<-colSums(y*w);yc<-sweep(y,2,ym,'-')
 beta<-solve(crossprod(x,x*w)+diag(lambda,ncol(x)),crossprod(x,yc*w))
 sweep(xx%*%beta,2,ym,'+')
}
weights<-function(z)1/as.numeric(table(z)[z])
select_predict<-function(train,test,features,target){
 ds<-unique(train$dataset);scores<-sapply(lambdas,function(lambda){
  mean(vapply(ds,function(s){a<-train[dataset!=s];b<-train[dataset==s]
   p<-ridge(a[,..features],a[,..target],b[,..features],weights(a$dataset),lambda)
   mean((p-as.matrix(b[,..target]))^2)
  },numeric(1)))
 })
 lambda<-lambdas[which.min(scores)]
 list(prediction=ridge(train[,..features],train[,..target],test[,..features],weights(train$dataset),lambda),lambda=lambda,inner_mse=min(scores))
}
predictions<-list();tuning<-list();metrics<-list();axis_metrics<-list();first_axis_predictions<-list()
for(rep in 1:10){
 d<-fread(file.path(out,paste0('analysis_units_partition',rep,'.tsv.gz')))
 d[,`:=`(log_n_treated=log(n_treated_A+n_treated_B),log_n_control=log(n_control_A+n_control_B),log_analytic=log(analytic_sampling+1e-4),log_risk=log(disagreement+1e-4))]
 for(a in axes)d[,(paste0('log_C_var_',a)):=log(get(paste0('C_var_',a))+1e-4)]
 state<-paste0('C_mean_',axes);variance<-paste0('log_C_var_',axes)
 support<-c('log_n_treated','log_n_control','C_mean_log_library')
 pilot<-paste0('response_',axes,'_A');targets<-paste0('response_',axes,'_B')
 feature_sets<-list(control_state=state,sampling_support=support,control_state_support=c(state,support),control_state_variability=c(state,variance,support),pilot_response=c(pilot,support),pilot_plus_state=c(pilot,state,variance,support))
 risk_sets<-list(sampling_support=support,control_state_support=c(state,support),control_state_variability=c(state,variance,support),pilot_response=c(pilot,support),sampling_aware=c('log_analytic',pilot,support),sampling_plus_state=c('log_analytic',pilot,state,support))
 for(s in unique(d$dataset)){
  train<-d[dataset!=s];test<-d[dataset==s];w<-weights(train$dataset);w<-w/sum(w)
  responses<-list(zero=matrix(0,nrow(test),length(axes)),training_mean=matrix(rep(colSums(as.matrix(train[,..targets])*w),each=nrow(test)),ncol=length(axes)),pilot_identity=as.matrix(test[,..pilot]))
  risks<-list(training_mean=rep(sum(train$disagreement*w),nrow(test)),analytic_sampling=test$analytic_sampling)
  for(m in names(feature_sets)){
   fit<-select_predict(train,test,feature_sets[[m]],targets);responses[[m]]<-fit$prediction
   tuning[[length(tuning)+1L]]<-data.table(partition_id=rep,held_out_study=s,task='response',model=m,lambda=fit$lambda,inner_mse=fit$inner_mse)
  }
  for(m in names(risk_sets)){
   fit<-select_predict(train,test,risk_sets[[m]],'log_risk');risks[[m]]<-pmax(exp(fit$prediction[,1])-1e-4,0)
   tuning[[length(tuning)+1L]]<-data.table(partition_id=rep,held_out_study=s,task='reliability',model=m,lambda=fit$lambda,inner_mse=fit$inner_mse)
  }
  for(m in names(responses)){
   p<-responses[[m]];err<-rowMeans((p-as.matrix(test[,..targets]))^2)
   ae<-(p-as.matrix(test[,..targets]))^2
   for(j in seq_along(axes))axis_metrics[[length(axis_metrics)+1L]]<-data.table(partition_id=rep,dataset=s,model=m,axis=axes[j],rmse=sqrt(mean(ae[,j])),n=nrow(test))
   if(rep==1L)for(j in seq_along(axes))first_axis_predictions[[length(first_axis_predictions)+1L]]<-cbind(test[,.(dataset,unit_id,context,perturbation)],data.table(model=m,axis=axes[j],observed=test[[targets[j]]],prediction=p[,j]))
   predictions[[length(predictions)+1L]]<-cbind(test[,.(partition_id,dataset,unit_id,context,biological_context,perturbation)],data.table(task='response',model=m,observed=test$disagreement,predicted=NA_real_,squared_error=err,log_squared_error=NA_real_))
   metrics[[length(metrics)+1L]]<-data.table(partition_id=rep,dataset=s,task='response',model=m,n=nrow(test),rmse=sqrt(mean(err)),spearman=NA_real_,log_rmse=NA_real_)
  }
  for(m in names(risks)){
   p<-risks[[m]];err<-(p-test$disagreement)^2;le<-(log(p+1e-4)-test$log_risk)^2
   predictions[[length(predictions)+1L]]<-cbind(test[,.(partition_id,dataset,unit_id,context,biological_context,perturbation)],data.table(task='reliability',model=m,observed=test$disagreement,predicted=p,squared_error=err,log_squared_error=le))
   metrics[[length(metrics)+1L]]<-data.table(partition_id=rep,dataset=s,task='reliability',model=m,n=nrow(test),rmse=sqrt(mean(err)),spearman=cor(p,test$disagreement,method='spearman'),log_rmse=sqrt(mean(le)))
  }
 }
 cat('Transfer partition',rep,'complete\n');flush.console()
}
p<-rbindlist(predictions);m<-rbindlist(metrics);t<-rbindlist(tuning)
fwrite(p,file.path(out,'leave_study_out_predictions.tsv.gz'),sep='\t')
fwrite(m,file.path(out,'leave_study_out_metrics.tsv'),sep='\t')
fwrite(t,file.path(out,'training_only_tuning.tsv'),sep='\t')
fwrite(rbindlist(axis_metrics),file.path(out,'response_metrics_by_axis.tsv'),sep='\t')
fwrite(rbindlist(first_axis_predictions),file.path(out,'response_predictions_by_axis_primary.tsv.gz'),sep='\t')
fwrite(m[,.(median_rmse=median(rmse),min_rmse=min(rmse),max_rmse=max(rmse),median_log_rmse=median(log_rmse),median_spearman=median(spearman)),by=.(dataset,task,model)],file.path(out,'partition_sensitivity.tsv'),sep='\t')
# Conditional intervention-cluster intervals for fixed first-partition predictions.
set.seed(202611062L);ci<-list();a<-p[partition_id==1L]
for(s in unique(a$dataset))for(task_name in c('response','reliability')){
 z<-a[dataset==s&task==task_name];ref<-if(task_name=='response')'pilot_identity'else'analytic_sampling'
 base<-z[model==ref,.(unit_id,reference_squared=squared_error,reference_log=log_squared_error)]
 for(me in setdiff(unique(z$model),ref)){
  v<-merge(z[model==me],base,by='unit_id');g<-split(seq_len(nrow(v)),v$perturbation)
  stat<-function(ii)if(task_name=='response')sqrt(mean(v$squared_error[ii]))-sqrt(mean(v$reference_squared[ii]))else sqrt(mean(v$log_squared_error[ii]))-sqrt(mean(v$reference_log[ii]))
  bo<-replicate(1000,stat(unlist(g[sample.int(length(g),length(g),replace=TRUE)],use.names=FALSE)))
  ci[[length(ci)+1L]]<-data.table(dataset=s,task=task_name,model=me,reference=ref,metric=if(task_name=='response')'response_RMSE_difference'else'log_disagreement_RMSE_difference',difference=stat(seq_len(nrow(v))),ci_low=quantile(bo,.025),ci_high=quantile(bo,.975),intervention_clusters=length(g),n_units=nrow(v),bootstrap_draws=length(bo),scope='fixed predictions and split; intervention clusters; not independent biological validation of baseline causality')
 }
}
fwrite(rbindlist(ci),file.path(out,'conditional_model_comparisons.tsv'),sep='\t')
