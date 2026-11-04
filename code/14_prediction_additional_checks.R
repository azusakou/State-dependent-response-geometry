#!/usr/bin/env Rscript
# Prediction validation from source-derived inputs.
suppressPackageStartupMessages({library(data.table);library(nnet);library(glmnet)})
run_nested_prediction_checks <- function(bundle, out, original_results=file.path(bundle,"results")) {
  bundle <- normalizePath(bundle)
  original_results <- normalizePath(original_results)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wr<-function(a,n)fwrite(a,file.path(out,n),sep='\t')
set.seed(19110313)
# Function-only loading cannot mutate primary code or result files.
fn_env<-new.env(parent=globalenv())
for(e in parse(file.path(bundle,'code/05_mixseq_models.R'))){
 if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,fn_env)
}
d<-as.data.table(readRDS(file.path(bundle,'inputs/mixseq_modeling.rds')))
r<-fread(file.path(original_results,'task_registry/state_transition_tasks.tsv'))
d<-d[pair_id%in%r[dataset_id=='mixseq_2020'&status=='eligible',task_id]];d[,row_id:=.I]
p<-fread(file.path(original_results,'mixseq_control_models/predictions.tsv'))
stopifnot(nrow(d)==120L,all(unique(p[,.(row_id,pair_id)])$pair_id==d$pair_id[unique(p[,.(row_id,pair_id)])$row_id]))
z_cols<-grep('^z[0-9]+_mean$',names(d),value=TRUE);y_cols<-paste0('p_',c('apoptosis','proliferation','stress'));fn_env$z_cols<-z_cols
methods<-c('local_response_ridge','kernel_transition','conditional_rate_network')
splits<-fread(file.path(original_results,'mixseq_control_models/split_assignments.tsv'))[split=='leave_cell_line']
fit_inner<-function(tr,te){
 stopifnot(!length(intersect(d$cell_line[tr],d$cell_line[te])))
 x1<-fn_env$design_matrix(d[tr],sort(unique(d$perturbation)));x2<-fn_env$design_matrix(d[te],sort(unique(d$perturbation)))
 xx<-fn_env$scale_train_test(x1,x2);zz<-fn_env$scale_train_test(as.matrix(d[tr,..z_cols]),as.matrix(d[te,..z_cols]));yy<-as.matrix(d[tr,..y_cols])
 unseen<-!d$perturbation[te]%in%d$perturbation[tr]
 cut<-fn_env$nearest_support_threshold(zz$train)
 near<-vapply(seq_along(te),function(i)min(sqrt(rowSums((zz$train-rep(zz$test[i,],each=length(tr)))^2))),numeric(1))
 rr<-list(local_response_ridge=fn_env$local_linear_predict(zz$train,xx$train,yy,zz$test,xx$test),kernel_transition=fn_env$graph_predict(zz$train,yy,zz$test,d$perturbation[tr],d$perturbation[te]),conditional_rate_network=fn_env$neural_ode_predict(xx$train,yy,ifelse(is.finite(d$time_hours[tr]),d$time_hours[tr],24),xx$test,ifelse(is.finite(d$time_hours[te]),d$time_hours[te],24)))
 rb<-rbindlist(lapply(methods,function(me){
  v<-rr[[me]];pred<-pmin(pmax(v$pred,0),1);hard<-!unseen&near<=cut;if(me=='kernel_transition')hard<-hard&v$support_n>=10
  data.table(row_id=te,pair_id=d$pair_id[te],cell_line=d$cell_line[te],method=me,hard_support=hard,uncertainty=v$uncertainty,mse=rowMeans((pred-as.matrix(d[te,..y_cols]))^2))
 }))
 rb
}
inners<-list();memberships<-list();thresholds<-list();decisions<-list()
# Risk budgets are exploratory tolerances on transcriptional endpoint-probability RMSE.
# They are not clinical targets, calibrated confidence levels or original prespecifications.
budgets<-c(.10,.15,.20)
for(of in sort(unique(splits$fold))){
 te<-splits[fold==of,row_id];tr<-setdiff(d$row_id,te)
 stopifnot(!length(intersect(d$cell_line[tr],d$cell_line[te])))
 inner_fold<-fn_env$balanced_group_folds(d$cell_line[tr],4L)
 memberships[[length(memberships)+1L]]<-data.table(outer_fold=of,row_id=d$row_id,pair_id=d$pair_id,cell_line=d$cell_line,role=ifelse(d$row_id%in%te,'outer_test','outer_train'),inner_validation_fold=ifelse(d$row_id%in%tr,inner_fold[match(d$row_id,tr)],NA_integer_))
 folds<-list()
 for(inf in 1:4){iv<-tr[inner_fold==inf];it<-setdiff(tr,iv);stopifnot(!length(intersect(c(it,iv),te)));v<-fit_inner(it,iv);v[,inner_fold:=inf];folds[[length(folds)+1L]]<-v}
 inn<-rbindlist(folds);inn[,outer_fold:=of];inners[[length(inners)+1L]]<-inn
 for(me in methods){
  a<-inn[method==me&hard_support&is.finite(uncertainty)]
  taus<-sort(unique(c(0,a$uncertainty,Inf)))
  cand<-rbindlist(lapply(taus,function(t){b<-a[uncertainty<=t];data.table(threshold=t,n_accepted=nrow(b),n_lines=uniqueN(b$cell_line),inner_rmse=if(nrow(b))sqrt(mean(b$mse))else NA_real_)}))
  # Outer errors/prediction labels are not read until the threshold is frozen here.
  for(budget in budgets){
   eligible<-cand[n_accepted>=10&n_lines>=3&inner_rmse<=budget]
   chosen<-if(nrow(eligible))eligible[order(-n_accepted,-threshold)][1]else data.table(threshold=-Inf,n_accepted=0L,n_lines=0L,inner_rmse=NA_real_)
   thresholds[[length(thresholds)+1L]]<-cbind(data.table(outer_fold=of,method=me,risk_budget=budget,inner_total=nrow(inn[method==me]),inner_hard_supported=nrow(a)),chosen)
   # Applicability is reconstructed from the outer training predictors only.
   z<-fn_env$scale_train_test(as.matrix(d[tr,..z_cols]),as.matrix(d[te,..z_cols]));distcut<-fn_env$nearest_support_threshold(z$train)
   nearest<-vapply(seq_along(te),function(i)min(sqrt(rowSums((z$train-rep(z$test[i,],each=length(tr)))^2))),numeric(1));np<-vapply(te,function(i)sum(d$perturbation[tr]==d$perturbation[i]),integer(1))
   outer<-p[split=='leave_cell_line'&fold==of&method==me][match(te,row_id)]
   hard<-np>0&nearest<=distcut;if(me=='kernel_transition')hard<-hard&np>=10
   outer[,`:=`(hard_support=hard,selected_uncertainty_threshold=chosen$threshold,risk_budget=budget,accepted=hard&is.finite(uncertainty)&uncertainty<=chosen$threshold)]
   # Check original decision reconstruction without changing it.
   stopifnot(all((hard&is.finite(outer$uncertainty)&outer$uncertainty<=.15)==!outer$abstain))
   for(e in c('apoptosis','proliferation','stress'))outer[,(paste0('se_',e)):=(get(paste0('predicted_',e))-get(paste0('observed_',e)))^2]
   outer[,mse:=rowMeans(.SD),.SDcols=paste0('se_',c('apoptosis','proliferation','stress'))]
   el<-p[split=='leave_cell_line'&fold==of&method=='elastic_net'][match(te,row_id)]
   outer[,elastic_mse:=rowMeans((as.matrix(el[,paste0('predicted_',c('apoptosis','proliferation','stress')),with=FALSE])-as.matrix(el[,paste0('observed_',c('apoptosis','proliferation','stress')),with=FALSE]))^2)]
   decisions[[length(decisions)+1L]]<-outer
  }
 }
 cat('Completed outer fold',of,'\n')
}
inn<-rbindlist(inners);mem<-rbindlist(memberships);sel<-rbindlist(thresholds);test<-rbindlist(decisions)
stopifnot(!anyDuplicated(test[,.(row_id,method,risk_budget)]),nrow(test)==120L*length(methods)*length(budgets))
wr(mem,'nested_calibration_membership.tsv');wr(inn,'nested_calibration_inner_predictions.tsv');wr(sel,'nested_selected_thresholds.tsv');wr(test,'nested_calibration_outer_predictions.tsv')
# Fixed-decision intervals on untouched outer predictions; NOT refitted nested-pipeline CIs.
B<-2000L;set.seed(19110314)
cluster_ci<-function(a){
 gs<-split(seq_len(nrow(a)),a$cell_line);ng<-length(gs)
 f<-function(x)c(response_rmse=sqrt(mean(x$mse)),elastic_rmse=sqrt(mean(x$elastic_mse)),response_minus_elastic=sqrt(mean(x$mse))-sqrt(mean(x$elastic_mse)))
 if(!nrow(a))return(data.table(metric=c('response_rmse','elastic_rmse','response_minus_elastic'),estimate=NA_real_,ci_low=NA_real_,ci_high=NA_real_))
 pt<-f(a);bo<-if(ng>1)replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);f(a[ii])})else matrix(NA_real_,3,B)
 rbindlist(lapply(seq_along(pt),function(j)data.table(metric=names(pt)[j],estimate=pt[j],ci_low=quantile(bo[j,],.025,na.rm=TRUE),ci_high=quantile(bo[j,],.975,na.rm=TRUE))))
}
summaries<-list()
for(me in methods)for(budget in budgets){
 a<-test[method==me&risk_budget==budget];b<-a[accepted==TRUE];ci<-cluster_ci(b)
 summaries[[length(summaries)+1L]]<-cbind(data.table(method=me,risk_budget=budget,n_total=nrow(a),n_accepted=nrow(b),n_lines=uniqueN(b$cell_line),coverage=nrow(b)/nrow(a),n_outer_folds_with_threshold=sel[method==me&risk_budget==budget,sum(threshold>-Inf)]),ci)
}
wr(rbindlist(summaries),'nested_calibration_summary.tsv')
wr(test[,.(n_total=.N,n_accepted=sum(accepted),rmse=if(any(accepted))sqrt(mean(mse[accepted]))else NA_real_,elastic_rmse=if(any(accepted))sqrt(mean(elastic_mse[accepted]))else NA_real_),by=.(method,risk_budget,fold)],'nested_calibration_by_outer_fold.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Nested threshold selection finished. Every outer case scored once per method/budget.\n')

}

run_norman_score_sensitivity <- function(bundle, out, original_results=file.path(bundle,"results")) {
  bundle <- normalizePath(bundle)
  original_results <- normalizePath(original_results)
  results <- original_results
dir.create(out,recursive=TRUE,showWarnings=FALSE);wr<-function(x,n)fwrite(x,file.path(out,n),sep='\t')
set.seed(19110315);B<-2000L
p<-fread(file.path(results,'grouped_validation/norman_predictions.tsv'));b<-fread(file.path(results,'extended_baselines/norman_baseline_predictions.tsv'))
p<-rbindlist(list(p[,.(split,fold,held_out,combination,method,predicted,observed,observed_class)],b[,.(split,fold,held_out,combination,method,predicted,observed,observed_class)]))
proj<-fread(file.path(bundle,'inputs/norman_alternative_projections.tsv'));meta<-fread(file.path(bundle,'inputs/norman_alternative_metadata.tsv'))
primary<-fread(file.path(results,'norman_control/combination_interactions.tsv'))
labels<-proj[,.(mean=mean(projection),se=sd(projection)/sqrt(.N),p=2*pt(-abs(mean(projection)/(sd(projection)/sqrt(.N))),.N-1L)),by=.(normalization,combination)]
labels[,q:=p.adjust(p,'BH'),by=normalization];labels<-merge(labels,meta,by=c('normalization','combination'))
labels[,positive:=q<.05&mean>0]
stopifnot(all(merge(labels[normalization=='log1p_CP10k'],primary[,.(combination,original_positive=interaction_class=='synergistic')],by='combination')[,positive==original_positive]))
# Recompute alternative exact sign-flip class without biological replication claim.
sf<-proj[normalization=='log1p_CP10k',{
 s<-as.matrix(expand.grid(rep(list(c(-1,1)),.N)));null<-as.numeric(s%*%projection)/.N
 .(mean=mean(projection),p=mean(abs(null)>=abs(mean(projection))-1e-12))
},by=combination];sf[,q:=p.adjust(p,'BH')];sf[,`:=`(normalization='CP10k_exact_technical_signflip',positive=q<.05&mean>0)]
labels<-rbindlist(list(labels[,.(normalization,combination,positive,mean,q)],sf[,.(normalization,combination,positive,mean,q)]))
wr(labels,'norman_label_sensitivity_inputs.tsv')
auc<-function(y,s){if(!any(y)||all(y))return(NA_real_);(sum(rank(s)[y])-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(!y))}
ap<-function(y,s){if(!any(y))return(NA_real_);o<-order(-s);yy<-y[o];ss<-s[o];i<-which(!duplicated(ss,fromLast=TRUE));sum((cumsum(yy)[i]/i)*diff(c(0,cumsum(yy)[i])))/sum(yy)}
metric<-function(a,score) c(auroc=auc(a$positive,a[[score]]),average_precision=ap(a$positive,a[[score]]),rmse=sqrt(mean((a[[score]]-a$observed)^2)),spearman=if(sd(a[[score]])>0&&sd(a$observed)>0)cor(a[[score]],a$observed,method='spearman')else NA_real_)
interval<-function(v){v<-v[is.finite(v)];if(length(v))unname(quantile(v,c(.025,.975)))else c(NA_real_,NA_real_)}
metrics<-list();paired<-list();methods<-c('nonlinear_xgboost','elastic_net','kernel_regression','differential_expression_ranking')
for(sp in unique(p$split))for(nm in unique(labels$normalization)){
 lab<-labels[normalization==nm];z<-merge(p[split==sp&method%in%methods],lab,by='combination')
 for(me in methods){a<-z[method==me];pt<-metric(a,'predicted')[1:2];gs<-split(seq_len(nrow(a)),a$held_out);ng<-length(gs)
 bo<-replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);metric(a[ii],'predicted')[1:2]})
 for(k in 1:2){ci<-interval(bo[k,]);metrics[[length(metrics)+1L]]<-data.table(split=sp,normalization=nm,method=me,n=nrow(a),positives=sum(a$positive),metric=names(pt)[k],estimate=pt[k],ci_low=ci[1],ci_high=ci[2],valid_bootstraps=sum(is.finite(bo[k,])),n_clusters=ng)}
 }
 for(me in setdiff(methods,'differential_expression_ranking')){
 a<-z[method==me];a<-merge(a,z[method=='differential_expression_ranking',.(combination,baseline_score=predicted)],by='combination');gs<-split(seq_len(nrow(a)),a$held_out);ng<-length(gs)
 pt<-metric(a,'predicted')-metric(a,'baseline_score')
 bo<-replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);metric(a[ii],'predicted')-metric(a[ii],'baseline_score')})
 # Continuous primary-response errors are unchanged under alternative class labels; report once only.
 ks<-if(nm=='log1p_CP10k')1:4 else 1:2
 for(k in ks){ci<-interval(bo[k,]);paired[[length(paired)+1L]]<-data.table(split=sp,normalization=nm,method=me,comparator='differential_expression_ranking',n=nrow(a),positives=sum(a$positive),metric=names(pt)[k],difference_model_minus_baseline=pt[k],ci_low=ci[1],ci_high=ci[2],valid_bootstraps=sum(is.finite(bo[k,])),n_clusters=ng)}
 }
}
wr(rbindlist(metrics),'norman_alternative_label_prediction_metrics.tsv');wr(rbindlist(paired),'norman_paired_baseline_differences.tsv')
cat('Completed paired baseline comparisons and frozen-score label sensitivity.\n')

}

run_raw_baseline_checks <- function(bundle='.',out='results/prediction_validation',original_results=file.path(bundle,'results')) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE);wr<-function(a,n)fwrite(a,file.path(out,n),sep='\t');set.seed(19110316)
 input_file<-file.path(bundle,'inputs/mixseq_raw_baseline_profiles.rds');if(!file.exists(input_file))input_file<-file.path(bundle,'mixseq_raw_baseline_profiles.rds')
 x<-readRDS(input_file);d<-as.data.table(x$modeling);d[,row_id:=.I];profiles<-x$expression
 stopifnot(nrow(d)==120,all(d$baseline_context%in%colnames(profiles)))
 s<-fread(file.path(original_results,'mixseq_control_models/split_assignments.tsv'))[split=='leave_cell_line']
 balanced<-function(g,k=5L){tab<-sort(table(g),decreasing=TRUE);load<-numeric(k);map<-integer(length(tab));names(map)<-names(tab);for(a in names(tab)){f<-which.min(load);map[a]<-f;load[f]<-load[f]+tab[[a]]};unname(map[as.character(g)])}
 scale_tt<-function(a,b){mu<-colMeans(a);ss<-apply(a,2,sd);ss[!is.finite(ss)|ss<1e-8]<-1;list(train=sweep(sweep(a,2,mu,'-'),2,ss,'/'),test=sweep(sweep(b,2,mu,'-'),2,ss,'/'))}
 ends<-c('apoptosis','proliferation','stress');ys<-paste0('p_',ends);all_predictions<-list();features<-list();membership<-list();diagnostics<-list()
 for(fo in sort(unique(s$fold))){
  te<-s[fold==fo,row_id];tr<-setdiff(d$row_id,te);stopifnot(!length(intersect(d$cell_line[tr],d$cell_line[te])))
  trctx<-sort(unique(d$baseline_context[tr]));tectx<-sort(unique(d$baseline_context[te]));stopifnot(!length(intersect(trctx,tectx)))
  # Feature ranking and PCA use unique outer-training baseline contexts only.
  g<-profiles[,trctx,drop=FALSE];v<-pmax(rowMeans(g*g)-rowMeans(g)^2,0);idx<-order(-v,rownames(g));idx<-head(idx[v[idx]>1e-12],2000L)
  ng<-length(idx);npcs<-min(30L,length(trctx)-1L,ng);fit<-prcomp(t(g[idx,,drop=FALSE]),center=TRUE,scale.=FALSE,rank.=npcs)
  zall<-predict(fit,t(profiles[idx,,drop=FALSE]))[,seq_len(npcs),drop=FALSE]
  z1<-zall[d$baseline_context[tr],,drop=FALSE];z2<-zall[d$baseline_context[te],,drop=FALSE];zz<-scale_tt(z1,z2)
  vocabulary<-sort(unique(d$perturbation));pert<-model.matrix(~factor(perturbation,levels=vocabulary)-1,data=d)
  tx<-log1p(ifelse(is.finite(d$time_hours),d$time_hours,0));xx<-scale_tt(cbind(z1,time=tx[tr],pert[tr,,drop=FALSE]),cbind(z2,time=tx[te],pert[te,,drop=FALSE]));y<-as.matrix(d[tr,..ys]);obs<-as.matrix(d[te,..ys])
  cvfold<-balanced(d$cell_line[tr],5L);enet<-matrix(NA_real_,length(te),3)
  for(j in 1:3){fitglm<-cv.glmnet(xx$train,y[,j],alpha=.5,foldid=cvfold,standardize=FALSE);enet[,j]<-as.numeric(predict(fitglm,xx$test,s='lambda.1se'))}
  # Same kernel and hard applicability definitions as the primary baseline-latent benchmark.
  dd<-as.matrix(dist(zz$train));diag(dd)<-Inf;supportcut<-as.numeric(quantile(apply(dd,1,min),.95));kernel<-matrix(NA_real_,length(te),3);unc<-numeric(length(te));np<-integer(length(te));near<-numeric(length(te))
  for(i in seq_along(te)){
   support<-which(d$perturbation[tr]==d$perturbation[te[i]]);np[i]<-length(support);if(!length(support))support<-seq_along(tr)
   near[i]<-min(sqrt(rowSums((zz$train-rep(zz$test[i,],each=length(tr)))^2)))
   d2<-rowSums((zz$train[support,,drop=FALSE]-rep(zz$test[i,],each=length(support)))^2);jj<-order(d2)[seq_len(min(60L,length(d2)))];ii<-support[jj];h<-sqrt(max(d2[jj],1e-8));w<-exp(-d2[jj]/(2*h*h));w<-w/sum(w);kernel[i,]<-colSums(y[ii,,drop=FALSE]*w);unc[i]<-sqrt(sum(w*rowSums((y[ii,,drop=FALSE]-rep(kernel[i,],each=length(ii)))^2))/3)
  }
  hard<-np>=10&near<=supportcut;accept<-hard&is.finite(unc)&unc<=.15
  predlist<-list(elastic_net=enet,kernel_transition=kernel)
  for(me in names(predlist)){
   pred<-pmin(pmax(predlist[[me]],0),1);p<-data.table(row_id=te,pair_id=d$pair_id[te],cell_line=d$cell_line[te],perturbation=d$perturbation[te],outer_fold=fo,method=me,accepted=if(me=='elastic_net')TRUE else accept,hard_support=hard,uncertainty=if(me=='elastic_net')NA_real_ else unc,mse=rowMeans((pred-obs)^2))
   for(j in 1:3){p[,(paste0('observed_',ends[j])):=obs[,j]];p[,(paste0('predicted_',ends[j])):=pred[,j]]}
   all_predictions[[length(all_predictions)+1L]]<-p
  }
  features[[length(features)+1L]]<-data.table(outer_fold=fo,gene=rownames(profiles)[idx],training_variance=v[idx])
  membership[[length(membership)+1L]]<-rbindlist(list(data.table(outer_fold=fo,baseline_context=trctx,role='fit_features_PCA'),data.table(outer_fold=fo,baseline_context=tectx,role='test_project_only')))
  diagnostics[[length(diagnostics)+1L]]<-data.table(outer_fold=fo,n_training_contexts=length(trctx),n_test_contexts=length(tectx),selected_genes=ng,n_PCs=npcs,training_distance_cutoff=supportcut,n_hard_supported=sum(hard),n_kernel_accepted=sum(accept),variance_explained=sum(fit$sdev[seq_len(npcs)]^2)/sum(fit$sdev^2))
 }
 p<-rbindlist(all_predictions);stopifnot(!anyDuplicated(p[,.(row_id,method)]),nrow(p)==240L)
 wr(p,'rawRNA_outer_predictions.tsv');wr(rbindlist(features),'rawRNA_training_selected_genes.tsv');wr(rbindlist(membership),'rawRNA_PCA_membership.tsv');wr(rbindlist(diagnostics),'rawRNA_fold_diagnostics.tsv')
 B<-2000L;set.seed(19110317);summary<-list()
 for(me in unique(p$method)){
  a<-p[method==me&accepted==TRUE];gs<-split(seq_len(nrow(a)),a$cell_line);ng<-length(gs)
  boot<-if(ng>1)replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);sqrt(mean(a$mse[ii]))})else NA_real_
  summary[[length(summary)+1L]]<-data.table(method=me,n_total=120,n_accepted=nrow(a),n_cell_lines=ng,coverage=nrow(a)/120,rmse=if(nrow(a))sqrt(mean(a$mse))else NA_real_,ci_low=if(ng>1)quantile(boot,.025)else NA_real_,ci_high=if(ng>1)quantile(boot,.975)else NA_real_)
 }
 wr(rbindlist(summary),'rawRNA_performance.tsv')
 # Both predictions evaluated on the exact kernel-accepted set.
 a<-merge(p[method=='kernel_transition'&accepted==TRUE,.(row_id,cell_line,kernel_mse=mse)],p[method=='elastic_net',.(row_id,elastic_mse=mse)],by='row_id')
 if(nrow(a)){
  gs<-split(seq_len(nrow(a)),a$cell_line);ng<-length(gs);boot<-if(ng>1)replicate(B,{ii<-unlist(gs[sample.int(ng,ng,replace=TRUE)],use.names=FALSE);sqrt(mean(a$kernel_mse[ii]))-sqrt(mean(a$elastic_mse[ii]))})else NA_real_
  wr(data.table(n=nrow(a),n_lines=ng,kernel_rmse=sqrt(mean(a$kernel_mse)),elastic_rmse=sqrt(mean(a$elastic_mse)),difference_kernel_minus_elastic=sqrt(mean(a$kernel_mse))-sqrt(mean(a$elastic_mse)),ci_low=if(ng>1)quantile(boot,.025)else NA_real_,ci_high=if(ng>1)quantile(boot,.975)else NA_real_),'rawRNA_matched_comparison.tsv')
 }else wr(data.table(n=0L,n_lines=0L,kernel_rmse=NA_real_,elastic_rmse=NA_real_,difference_kernel_minus_elastic=NA_real_,ci_low=NA_real_,ci_high=NA_real_),'rawRNA_matched_comparison.tsv')
 cat('Outer train-fold RNA selection/PCA benchmark finished.\n')
}

run_prediction_validation <- function(bundle=".", out=file.path(bundle,"results/prediction_validation"), original_results=file.path(bundle,"results")) {
 run_nested_prediction_checks(bundle,out,original_results)
 run_norman_score_sensitivity(bundle,out,original_results)
 run_raw_baseline_checks(bundle,out,original_results)
 run_growth_bridge(bundle,out,original_results)
 invisible(out)
}


# Orthogonal bulk-growth assay check; no inference of stable fates or independent pair replicates.
run_growth_bridge <- function(bundle='.',out=file.path(bundle,'results/prediction_validation'),original_results=file.path(bundle,'results')) {
 suppressPackageStartupMessages(library(data.table));dir.create(out,recursive=TRUE,showWarnings=FALSE)
 a<-fread(file.path(bundle,'inputs/norman_growth_fitness.tsv'));stopifnot(!anyDuplicated(a$combination),nrow(a)==131L)
 a[,matched:=is.finite(growth_rep1)&is.finite(growth_rep2)];a[,growth_average:=ifelse(matched,(growth_rep1+growth_rep2)/2,NA_real_)]
 t<-fread(file.path(original_results,'task_registry/state_transition_tasks.tsv'))[dataset_id=='norman_2019_combinatorial'];t[,combination:=sub('^norman_k562\\|','',task_id)]
 a<-merge(a,t[,.(combination,basin=target_state_id,task_status=status)],by='combination',all.x=TRUE);stopifnot(!anyNA(a$basin))
 z<-a[matched==TRUE&task_status=='eligible'];stopifnot(nrow(z)==110L,uniqueN(z$basin)==15L)
 b<-z[,.(n_pairs=.N,growth_median=median(growth_average),growth_q25=as.numeric(quantile(growth_average,.25)),growth_q75=as.numeric(quantile(growth_average,.75)),growth_min=min(growth_average),growth_max=max(growth_average),rep1_median=median(growth_rep1),rep2_median=median(growth_rep2)),by=basin][order(basin)]
 eta<-function(y,g)sum((ave(y,g,FUN=mean)-mean(y))^2)/sum((y-mean(y))^2)
 m<-rbindlist(lapply(c('growth_average','growth_rep1','growth_rep2'),function(nm){y<-z[[nm]];dd<-abs(outer(y,y,'-'));same<-outer(z$basin,z$basin,'==');up<-upper.tri(dd);data.table(measure=nm,n_pairs=nrow(z),n_basins=uniqueN(z$basin),n_basins_single_pair=sum(table(z$basin)==1),eta_squared=eta(y,z$basin),within_comparisons=sum(same&up),between_comparisons=sum(!same&up),median_within_distance=median(dd[same&up]),median_between_distance=median(dd[!same&up]))}))
 support<-data.table(all_pairs=nrow(a),matched_pairs=sum(a$matched),eligible_pairs=sum(a$task_status=='eligible'),matched_eligible_pairs=nrow(z),eligible_basins=uniqueN(a[task_status=='eligible',basin]),matched_eligible_basins=uniqueN(z$basin),replicate_pearson=cor(z$growth_rep1,z$growth_rep2),replicate_spearman=cor(z$growth_rep1,z$growth_rep2,method='spearman'))
 fwrite(a,file.path(out,'functional_condition_mapping.tsv'),sep='\t');fwrite(b,file.path(out,'basin_growth_summary.tsv'),sep='\t');fwrite(m,file.path(out,'functional_descriptive_metrics.tsv'),sep='\t');fwrite(support,file.path(out,'functional_match_summary.tsv'),sep='\t')
 invisible(support)
}

run_prediction_validation()
