#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
out<-'paper19opp-revision2/reconstruction/independent_lineage'
started<-Sys.time();set.seed(20261106)
cache<-readRDS(file.path(out,'analysis_units.rds'));features<-cache$features
lambdas<-10^seq(-4,4,length.out=17)
mean_y<-function(y,w)colSums(y*w)/sum(w)
sd_y<-function(y,w){mu<-mean_y(y,w);pmax(sqrt(colSums(sweep(y,2,mu,'-')^2*w)/sum(w)),1e-8)}
ridge_grid<-function(x,y,newx,w,lambda=lambdas){
 keep<-w>0;x<-x[keep,,drop=FALSE];y<-y[keep,,drop=FALSE];w<-w[keep];w<-w/sum(w)
 xm<-colSums(x*w);xx<-sweep(x,2,xm,'-');xs<-sqrt(colSums(xx^2*w));xs[!is.finite(xs)|xs<1e-8]<-1
 xx<-sweep(xx,2,xs,'/');tx<-sweep(sweep(newx,2,xm,'-'),2,xs,'/')
 ym<-colSums(y*w);yy<-sweep(y,2,ym,'-')
 ee<-eigen(crossprod(xx*sqrt(w)),symmetric=TRUE);ev<-pmax(ee$values,0)
 rot<-crossprod(ee$vectors,crossprod(xx*w,yy));proj<-tx%*%ee$vectors
 lapply(lambda,function(ll)sweep(proj%*%sweep(rot,1,ev+ll,'/'),2,ym,'+'))
}
condition_prediction<-function(y,condition,w,testcondition,parent_train=NULL,parent_test=NULL){
 if(!is.null(parent_train))y<-y-parent_train
 mu<-mean_y(y,w);pred<-matrix(rep(mu,each=length(testcondition)),nrow=length(testcondition))
 for(co in unique(condition)){ii<-which(condition==co&w>0);jj<-which(testcondition==co);if(length(ii)&&length(jj))pred[jj,]<-matrix(rep(mean_y(y[ii,,drop=FALSE],w[ii]),each=length(jj)),nrow=length(jj))}
 if(!is.null(parent_test))pred<-pred+parent_test
 pred
}
risk_numerator<-function(pred,y,w,yscale)sum(w*rowMeans(sweep(pred-y,2,yscale,'/')^2))
run_cv<-function(d,xb,xs,y,yp,w,outer,inner,save_predictions=FALSE){
 models<-c('condition_mean','counts_ridge','parent_persistence','selected_simple','state_ridge');pred<-lapply(models,function(z)matrix(NA_real_,nrow(d),3));names(pred)<-models
 scales<-matrix(NA_real_,nrow(d),3);decisions<-list()
 for(k in 1:5){
  test<-which(outer==k);train<-which(outer!=k);if(sum(w[test])<=0)next
  bc<-numeric(length(lambdas));sc<-bc;simple<-c(condition_mean=0,parent_persistence=0);total<-0
  for(j in 1:4){
   va<-train[inner[[k]][train]==j];tr<-train[inner[[k]][train]!=j];if(sum(w[va])<=0||sum(w[tr])<=0)next
   scale<-sd_y(y[tr,,drop=FALSE],w[tr]);total<-total+sum(w[va])
   bp<-ridge_grid(xb[tr,,drop=FALSE],y[tr,,drop=FALSE],xb[va,,drop=FALSE],w[tr]);sp<-ridge_grid(xs[tr,,drop=FALSE],y[tr,,drop=FALSE],xs[va,,drop=FALSE],w[tr])
   for(h in seq_along(lambdas)){bc[h]<-bc[h]+risk_numerator(bp[[h]],y[va,,drop=FALSE],w[va],scale);sc[h]<-sc[h]+risk_numerator(sp[[h]],y[va,,drop=FALSE],w[va],scale)}
   mm<-condition_prediction(y[tr,,drop=FALSE],d$sequence[tr],w[tr],d$sequence[va]);pp<-condition_prediction(y[tr,,drop=FALSE],d$sequence[tr],w[tr],d$sequence[va],yp[tr,,drop=FALSE],yp[va,,drop=FALSE])
   simple['condition_mean']<-simple['condition_mean']+risk_numerator(mm,y[va,,drop=FALSE],w[va],scale)
   simple['parent_persistence']<-simple['parent_persistence']+risk_numerator(pp,y[va,,drop=FALSE],w[va],scale)
  }
  stopifnot(total>0)
  ib<-which.min(bc);is<-which.min(sc);simple<-c(simple,counts_ridge=bc[ib]);winner<-names(which.min(simple))
  pred$condition_mean[test,]<-condition_prediction(y[train,,drop=FALSE],d$sequence[train],w[train],d$sequence[test])
  pred$parent_persistence[test,]<-condition_prediction(y[train,,drop=FALSE],d$sequence[train],w[train],d$sequence[test],yp[train,,drop=FALSE],yp[test,,drop=FALSE])
  pred$counts_ridge[test,]<-ridge_grid(xb[train,,drop=FALSE],y[train,,drop=FALSE],xb[test,,drop=FALSE],w[train],lambdas[ib])[[1]]
  pred$state_ridge[test,]<-ridge_grid(xs[train,,drop=FALSE],y[train,,drop=FALSE],xs[test,,drop=FALSE],w[train],lambdas[is])[[1]]
  pred$selected_simple[test,]<-pred[[winner]][test,,drop=FALSE]
  scales[test,]<-matrix(rep(sd_y(y[train,,drop=FALSE],w[train]),each=length(test)),nrow=length(test))
  decisions[[k]]<-data.table(fold=k,selected_simple=winner,counts_lambda=lambdas[ib],state_lambda=lambdas[is],counts_inner_rmse=sqrt(bc[ib]/total),state_inner_rmse=sqrt(sc[is]/total),simple_inner_rmse=sqrt(min(simple)/total))
 }
 active<-which(w>0);stopifnot(all(is.finite(scales[active,,drop=FALSE])))
 metrics<-rbindlist(lapply(models,function(m){pr<-pred[[m]];stopifnot(all(is.finite(pr[active,,drop=FALSE])));error<-pr-y
 data.table(model=m,scaled_rmse=sqrt(sum(w[active]*rowMeans((error[active,,drop=FALSE]/scales[active,,drop=FALSE])^2))/sum(w)),rmse_apoptosis=sqrt(sum(w[active]*error[active,1]^2)/sum(w)),rmse_proliferation=sqrt(sum(w[active]*error[active,2]^2)/sum(w)),rmse_stress=sqrt(sum(w[active]*error[active,3]^2)/sum(w)))}))
 list(metrics=metrics,predictions=if(save_predictions)pred else NULL,scales=if(save_predictions)scales else NULL,decisions=if(save_predictions)rbindlist(decisions)else NULL)
}
allmetrics<-list();bootstrap<-list()
for(setting in c('primary','strict')){
 gate<-if(setting=='primary')c(3,2)else c(10,5)
 d<-copy(cache$units[parent_n_cells>=gate[1]&desc_n_cells>=gate[2]])
 ids<-sort(unique(d$lineage));set.seed(20261106+match(setting,c('primary','strict')))
 fmap<-setNames(sample(rep(1:5,length.out=length(ids))),ids);outer<-unname(fmap[d$lineage]);inner<-vector('list',5)
 for(k in 1:5){trainids<-sort(unique(d$lineage[outer!=k]));set.seed(20261106+100*k+match(setting,c('primary','strict')));im<-setNames(sample(rep(1:4,length.out=length(trainids))),trainids);inner[[k]]<-unname(im[d$lineage])}
 d[,outer_fold:=outer];d[,weight:=1/.N,by=lineage]
 y<-as.matrix(d[,paste0('desc_score_',c('apoptosis','proliferation','stress')),with=FALSE]);yp<-as.matrix(d[,paste0('parent_score_',c('apoptosis','proliferation','stress')),with=FALSE])
 factors<-model.matrix(~0+factor(sequence),d);counts<-cbind(log_parent=log1p(d$parent_n_cells),log_desc=log1p(d$desc_n_cells),parent_log_library=d$parent_mean_log_library,parent_detected_genes=d$parent_mean_detected_genes)
 xb<-cbind(factors,counts);xs<-cbind(xb,as.matrix(d[,paste0('parent_',features),with=FALSE]));stopifnot(all(is.finite(c(xb,xs,y,yp))))
 fwrite(d[,.(lineage,first,second,parent_n_cells,desc_n_cells,outer_fold,weight)],file.path(out,paste0(setting,'_folds.tsv')),sep='\t')
 fit<-run_cv(d,xb,xs,y,yp,d$weight,outer,inner,TRUE);fit$metrics[,`:=`(setting=setting,tasks=nrow(d),lineages=length(ids))];allmetrics[[setting]]<-fit$metrics
 fwrite(fit$decisions,file.path(out,paste0(setting,'_training_decisions.tsv')),sep='\t')
 predictions<-rbindlist(lapply(names(fit$predictions),function(m){z<-copy(d[,.(lineage,first,second,outer_fold,weight)]);z[,model:=m];for(j in 1:3){nm<-c('apoptosis','proliferation','stress')[j];z[,(paste0('observed_',nm)):=y[,j]];z[,(paste0('predicted_',nm)):=fit$predictions[[m]][,j]];z[,(paste0('training_sd_',nm)):=fit$scales[,j]]};z}))
 fwrite(predictions,file.path(out,paste0(setting,'_heldout_predictions.tsv.gz')),sep='\t')
 cat(setting,'point estimates, tasks',nrow(d),'lineages',length(ids),'\n');print(fit$metrics);flush.console()
 set.seed(20261106+10000+match(setting,c('primary','strict')))
 for(b in 1:249){
  sampled<-table(factor(sample(ids,length(ids),replace=TRUE),levels=ids));w<-d$weight*as.numeric(sampled[match(d$lineage,ids)])
  ff<-run_cv(d,xb,xs,y,yp,w,outer,inner,FALSE)$metrics;ff[,`:=`(setting=setting,replicate=b)];bootstrap[[paste(setting,b)]]<-ff
  if(b%%25==0){fwrite(rbindlist(bootstrap),file.path(out,'bootstrap_metrics_checkpoint.tsv'),sep='\t');cat(setting,'bootstrap',b,'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'seconds\n');flush.console()}
 }
}
metrics<-rbindlist(allmetrics);bs<-rbindlist(bootstrap);fwrite(metrics,file.path(out,'heldout_metrics.tsv'),sep='\t');fwrite(bs,file.path(out,'bootstrap_metrics.tsv'),sep='\t')
comparisons<-rbindlist(lapply(c('primary','strict'),function(s){rbindlist(lapply(c('selected_simple','condition_mean','counts_ridge','parent_persistence'),function(m){a<-bs[setting==s&model=='state_ridge',.(replicate,state=scaled_rmse)];b<-bs[setting==s&model==m,.(replicate,simple=scaled_rmse)];zz<-merge(a,b,by='replicate');diff<-zz$state-zz$simple;point<-metrics[setting==s&model=='state_ridge',scaled_rmse]-metrics[setting==s&model==m,scaled_rmse];data.table(setting=s,comparison=paste0('state_ridge_minus_',m),difference=point,lower=quantile(diff,.025),upper=quantile(diff,.975),bootstrap_draws=length(diff),negative_fraction=mean(diff<0))}))}))
fwrite(comparisons,file.path(out,'paired_model_comparisons.tsv'),sep='\t');print(comparisons)
writeLines(c(paste('elapsed_seconds',as.numeric(difftime(Sys.time(),started,units='secs'))),capture.output(sessionInfo())),file.path(out,'prediction_session_info.txt'))
cat('Completed',format(Sys.time(),tz='UTC'),'UTC\n')
