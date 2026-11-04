#!/usr/bin/env Rscript
# Strict nested downstream prediction using separate train-only scVI fits.
# Usage: Rscript nested_scvi_models.R BUNDLE_ROOT LATENT_ROOT OUTPUT_DIR [--primary-only] [--nn-maxit=600]
# Default also scores known-vocabulary low-support conditions with frozen primary policies.
# This stage never substitutes the original globally fitted latent coordinates.
suppressPackageStartupMessages({library(data.table);library(glmnet);library(nnet)})
run_nested_scvi_models<-function(bundle,latent_root,out,include_low_support=TRUE,nn_maxit=600L){
stopifnot(length(nn_maxit)==1L,is.finite(nn_maxit),nn_maxit>0,nn_maxit==as.integer(nn_maxit));nn_maxit<-as.integer(nn_maxit)
bundle<-normalizePath(bundle);latent_root<-normalizePath(latent_root)
latent_table<-NULL;context_table<-NULL
if(file.exists(latent_root)&&!dir.exists(latent_root)){
 cache_latent<-readRDS(latent_root)
 if(is.list(cache_latent)&&'context_latents'%in%names(cache_latent))context_table<-as.data.table(cache_latent$context_latents)
 if(is.list(cache_latent)&&'latents'%in%names(cache_latent))cache_latent<-cache_latent$latents
 if(is.data.frame(cache_latent))latent_table<-as.data.table(cache_latent)
 else if(is.list(cache_latent))latent_table<-rbindlist(cache_latent,use.names=TRUE,fill=TRUE)
 else stop('Consolidated latent RDS must contain a table or list of tables.')
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,name)fwrite(x,file.path(out,name),sep='\t')
fn<-new.env(parent=globalenv())
for(e in parse(file.path(bundle,'code','05_mixseq_models.R'))){
 if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,fn)
}
# Audit the unchanged solvers without additional fits or RNG calls.  The copied
# original NN function adds only a recorder and an explicit iteration-budget parameter;
# default 600 reproduces the original call, while 6000 is a prespecified sensitivity.
optimizer_rows<-list();fn$nn_iteration_limit<-nn_maxit
fn$record_nnet_optimizer<-function(fit,seed,n_training,n_query){
 cx<-fn$fit_audit_context
 optimizer_rows[[length(optimizer_rows)+1L]]<<-data.table(fit_id=cx$fit_id,outer_fold=cx$outer_fold,inner_fold=cx$inner_fold,method='conditional_rate_network',endpoint=NA_character_,seed=as.integer(seed),n_training=as.integer(n_training),n_query=as.integer(n_query),iteration_limit=fn$nn_iteration_limit,reported_convergence=as.integer(fit$convergence),objective=as.numeric(fit$value),objective_finite=is.finite(fit$value),npasses=NA_integer_,jerr=NA_integer_,requested_lambda_count=NA_integer_,fitted_lambda_count=NA_integer_,requested_lambda_min=NA_real_,requested_lambda_max=NA_real_,fitted_lambda_min=NA_real_,fitted_lambda_max=NA_real_,lambda_grid_complete=NA,parameters_finite=all(is.finite(fit$wts)),status=if(fit$convergence==0L)'converged'else if(fit$convergence==1L)'iteration_limit_reached'else 'other_reported_status')
 invisible(NULL)
}
nnet_audit_insertions<-0L
instrument_nnet_body<-function(expr){
 if(!is.call(expr))return(expr)
 if(identical(expr[[1L]],as.name('{'))){
  parts<-as.list(expr)[-1L];updated<-list(as.name('{'))
  for(part in parts){
   is_nnet_fit<-is.call(part)&&identical(part[[1L]],as.name('<-'))&&identical(part[[2L]],as.name('fit'))&&is.call(part[[3L]])&&identical(part[[3L]][[1L]],quote(nnet::nnet))
   if(is_nnet_fit){
    stopifnot(identical(as.numeric(part[[3L]][['maxit']]),600))
    part[[3L]][['maxit']]<-quote(nn_iteration_limit)
   }
   updated[[length(updated)+1L]]<-instrument_nnet_body(part)
   if(is_nnet_fit){
    updated[[length(updated)+1L]]<-quote(record_nnet_optimizer(fit,1904+seeds[s],nrow(x_train),nrow(x_test)))
    nnet_audit_insertions<<-nnet_audit_insertions+1L
   }
  }
  return(as.call(updated))
 }
 as.call(lapply(as.list(expr),instrument_nnet_body))
}
body(fn$neural_ode_predict)<-instrument_nnet_body(body(fn$neural_ode_predict))
stopifnot(nnet_audit_insertions==1L)

all_tasks<-as.data.table(readRDS(file.path(bundle,'inputs','mixseq_modeling.rds')))
all_tasks[,baseline_context:=paste(experiment,ifelse(is.na(time_hours),'NA',as.character(time_hours)),cell_line,sep='|')]
registry<-fread(file.path(bundle,'tables','task_catalogue.tsv.gz'))
tasks<-copy(all_tasks[pair_id%in%registry[dataset_id=='mixseq_2020'&status=='eligible',task_id]])
tasks[,row_id:=.I]
tasks[,baseline_context:=paste(experiment,ifelse(is.na(time_hours),'NA',as.character(time_hours)),cell_line,sep='|')]
ends<-c('apoptosis','proliferation','stress');y_cols<-paste0('p_',ends)
z_cols<-paste0('z',0:29,'_mean');fn$z_cols<-z_cols
stopifnot(nrow(tasks)==120L,!anyDuplicated(tasks$pair_id),all(complete.cases(tasks[,..y_cols])))
# Intervention names are fixed experimental metadata, not outcome-derived features.
perturbation_vocabulary<-sort(unique(tasks$perturbation))
splits<-fread(file.path(bundle,'tables','mixseq_validation_splits.tsv.gz'))[split=='leave_cell_line']
stopifnot(nrow(splits)==120L,setequal(splits$row_id,tasks$row_id))
# The selection sensitivity never contributes training labels, vocabulary, lambda or cutoffs.
# All extra conditions are on benchmark lines; novel identities are audited, not encoded.
extra<-copy(all_tasks[cell_line%in%tasks$cell_line&!pair_id%in%tasks$pair_id])
extra[,row_id:=120L+.I]
line_folds<-unique(cbind(tasks[,.(cell_line)],data.table(outer_fold=splits$fold[match(tasks$row_id,splits$row_id)])))
stopifnot(!anyDuplicated(line_folds$cell_line))
extra[,outer_fold:=line_folds$outer_fold[match(cell_line,line_folds$cell_line)]]
extra[,`:=`(known_intervention=perturbation%in%perturbation_vocabulary,
 minimum_support=pmin(baseline_cells,treated_cells),
 maximum_program_se=pmax(p_apoptosis_se,p_proliferation_se,p_stress_se))]
extra[,`:=`(count_bin=ifelse(minimum_support<50L,'20-49','50-99'),
 se_bin=ifelse(maximum_program_se<=.05,'<=0.05',ifelse(maximum_program_se<=.10,'(0.05,0.10]','>0.10')),
 support_deficit=ifelse(baseline_cells<100L&treated_cells<100L,'both',ifelse(baseline_cells<100L,'baseline_only','treated_only')))]
stopifnot(nrow(extra)==210L,sum(extra$known_intervention)==144L,all(extra$minimum_support>=20L&extra$minimum_support<100L),all(complete.cases(extra[,..y_cols])),all(is.finite(extra$maximum_program_se)),!anyDuplicated(extra$pair_id),!anyNA(extra$outer_fold))
extra_scope<-extra[,.(row_id,pair_id,cell_line,baseline_context,perturbation,outer_fold,baseline_cells,treated_cells,minimum_support,maximum_program_se,count_bin,se_bin,support_deficit,known_intervention)]
extra_scope[,`:=`(numerically_evaluated=known_intervention,
 prespecified_response_identity_rejection=!known_intervention,
 evaluation_scope=ifelse(known_intervention,'extra low-support conditions; frozen primary-panel model/policy','outside fixed primary intervention vocabulary; no static or response prediction computed'))]

lambda_grid<-10^seq(0,-4,length.out=17L)
risk_budgets<-c(.10,.15,.20)
response_methods<-c('local_response_ridge','kernel_transition','conditional_rate_network')
all_methods<-c('training_mean','ridge_static','elastic_net',response_methods)
manifest_parts<-list()
for(of in sort(unique(splits$fold))){
 test<-splits[fold==of,row_id];train<-setdiff(tasks$row_id,test)
 stopifnot(!length(intersect(tasks$cell_line[train],tasks$cell_line[test])))
 inner<-fn$balanced_group_folds(tasks$cell_line[train],4L)
 m<-tasks[,.(row_id,pair_id,cell_line,baseline_context,perturbation,time_hours)]
 m[,`:=`(fit_id=paste0('outer_',of),outer_fold=of,inner_fold=0L,task_role=ifelse(row_id%in%train,'train','test'))]
 manifest_parts[[length(manifest_parts)+1L]]<-m
 for(inf in 1:4){
  valid<-train[inner==inf]
  m<-tasks[,.(row_id,pair_id,cell_line,baseline_context,perturbation,time_hours)]
  m[,`:=`(fit_id=paste0('outer_',of,'_inner_',inf),outer_fold=of,inner_fold=inf,task_role=ifelse(row_id%in%test,'outer_test_excluded',ifelse(row_id%in%valid,'validation','train')))]
  manifest_parts[[length(manifest_parts)+1L]]<-m
 }
}
manifest<-rbindlist(manifest_parts)
write_tsv(manifest,'nested_scvi_fold_manifest.tsv')
write_tsv(data.table(lambda=lambda_grid,alpha=.5,selection='minimum pooled inner out-of-fold endpoint MSE; ties prefer larger lambda'),'nested_scvi_lambda_grid.tsv')
write_tsv(data.table(nn_iteration_limit=nn_maxit,nn_seeds='1905,1906,1907',lambda_selection='grouped inner minimum MSE',threshold_selection='grouped inner empirical budgets 0.10,0.15,0.20',minimum_inner_accepted=10L,minimum_inner_accepted_lines=3L,low_support_queries=include_low_support),'nested_scvi_protocol.tsv')

fit_table_checks<-list()
read_latent_fit<-function(fit_id){
 fit_key<-fit_id
 if(is.null(latent_table)){
  f<-file.path(latent_root,'fits',fit_id,'task_latent_means.tsv')
  if(!file.exists(f))stop('Missing trained scVI fit output: ',f)
  z<-fread(f)
 }else z<-copy(latent_table[fit_id==fit_key])
 required<-c('fit_id','outer_fold','inner_fold','row_id','pair_id','cell_line','baseline_context','task_role',z_cols)
 stopifnot(all(required%in%names(z)),!anyDuplicated(z$pair_id),all(z$fit_id==fit_id),all(is.finite(as.matrix(z[,..z_cols]))))
 fit_key<-fit_id
 expected<-manifest[fit_id==fit_key&task_role!='outer_test_excluded']
 stopifnot(setequal(z$pair_id,expected$pair_id),!any(z$task_role=='outer_test_excluded'))
 z<-z[match(expected$pair_id,pair_id)]
 stopifnot(identical(z$pair_id,expected$pair_id),all(z$row_id==expected$row_id),all(z$task_role==expected$task_role),all(z$cell_line==expected$cell_line),all(z$baseline_context==expected$baseline_context))
 d<-copy(tasks[match(z$pair_id,pair_id)])
 d[,(z_cols):=z[,z_cols,with=FALSE]]
 d[,task_role:=z$task_role]
 stopifnot(!length(intersect(d[task_role=='train']$cell_line,d[task_role!='train']$cell_line)),!length(intersect(d[task_role=='train']$baseline_context,d[task_role!='train']$baseline_context)))
 fit_table_checks[[length(fit_table_checks)+1L]]<<-data.table(fit_id=fit_id,outer_fold=unique(z$outer_fold),inner_fold=unique(z$inner_fold),training_tasks=sum(d$task_role=='train'),query_tasks=sum(d$task_role!='train'),training_lines=uniqueN(d[task_role=='train']$cell_line),query_lines=uniqueN(d[task_role!='train']$cell_line),task_identifier_and_role_match=TRUE,line_and_context_disjoint=TRUE,no_outer_test_in_inner_output=TRUE)
 list(data=d,fit_id=fit_id,outer_fold=unique(z$outer_fold),inner_fold=unique(z$inner_fold))
}

read_extra_contexts<-function(fit_id,outer_fold){
 if(is.null(context_table)){
  if(!dir.exists(latent_root))stop('Consolidated latent RDS has no context_latents for requested low-support evaluation.')
  context_path<-file.path(latent_root,'fits',fit_id,'all_baseline_context_latent_means.tsv')
  if(!file.exists(context_path))stop('Missing all-context outer encoder output: ',context_path)
  cz<-fread(context_path)
 }else{fit_key<-fit_id;cz<-copy(context_table[fit_id==fit_key])}
 required<-c('baseline_context','cell_line','baseline_ncells',z_cols)
 stopifnot(all(required%in%names(cz)),!anyDuplicated(cz$baseline_context))
 fold_key<-outer_fold;ex<-copy(extra[outer_fold==fold_key])
 ix<-match(ex$baseline_context,cz$baseline_context)
 stopifnot(!anyNA(ix),all(cz$cell_line[ix]==ex$cell_line),all(cz$baseline_ncells[ix]==ex$baseline_cells),all(is.finite(as.matrix(cz[ix,..z_cols]))))
 ex[,(z_cols):=cz[ix,z_cols,with=FALSE]]
 ex[,task_role:='extra_test']
 # Every additional condition is encoded once by the encoder that excludes its line.
 ex[known_intervention==TRUE]
}

predict_fit<-function(fit_data,selected_lambda=NULL){
 d<-fit_data$data;tr<-which(d$task_role=='train');te<-which(d$task_role!='train')
 fn$fit_audit_context<-list(fit_id=fit_data$fit_id,outer_fold=fit_data$outer_fold,inner_fold=fit_data$inner_fold)
 stopifnot(length(tr)>=20L,length(te)>0L)
 xtr<-fn$design_matrix(d[tr],perturbation_vocabulary);xte<-fn$design_matrix(d[te],perturbation_vocabulary)
 xx<-fn$scale_train_test(xtr,xte);zz<-fn$scale_train_test(as.matrix(d[tr,..z_cols]),as.matrix(d[te,..z_cols]))
 stopifnot(all(is.finite(xx$train)),all(is.finite(xx$test)),all(is.finite(zz$train)),all(is.finite(zz$test)))
 yy<-as.matrix(d[tr,..y_cols]);ys<-as.matrix(d[te,..y_cols])
 unseen<-!d$perturbation[te]%in%d$perturbation[tr]
 dist_cut<-as.numeric(fn$nearest_support_threshold(zz$train))
 nearest<-vapply(seq_along(te),function(i)min(sqrt(rowSums((zz$train-rep(zz$test[i,],each=length(tr)))^2))),numeric(1))
 rr<-list()
 rr$training_mean<-list(pred=matrix(rep(colMeans(yy),each=length(te)),nrow=length(te)),uncertainty=rep(mean(apply(yy,2,sd)),length(te)))
 ridge<-fn$ridge_predict(xx$train,yy,xx$test,lambda=2)
 rr$ridge_static<-list(pred=ridge,uncertainty=rep(sqrt(mean((yy-fn$ridge_predict(xx$train,yy,xx$train,lambda=2))^2)),length(te)))
 rr$local_response_ridge<-fn$local_linear_predict(zz$train,xx$train,yy,zz$test,xx$test)
 rr$kernel_transition<-fn$graph_predict(zz$train,yy,zz$test,d$perturbation[tr],d$perturbation[te])
 rr$conditional_rate_network<-fn$neural_ode_predict(xx$train,yy,ifelse(is.finite(d$time_hours[tr]),d$time_hours[tr],24),xx$test,ifelse(is.finite(d$time_hours[te]),d$time_hours[te],24))
 enet_array<-array(NA_real_,dim=c(length(te),3L,length(lambda_grid)))
 enet_selected<-matrix(NA_real_,length(te),3L)
 for(j in 1:3){
  fit<-glmnet::glmnet(xx$train,yy[,j],alpha=.5,lambda=lambda_grid,standardize=FALSE,thresh=1e-10,maxit=100000L)
  full_grid<-length(fit$lambda)==length(lambda_grid)&&isTRUE(all.equal(as.numeric(fit$lambda),as.numeric(lambda_grid),tolerance=1e-12))
  optimizer_rows[[length(optimizer_rows)+1L]]<<-data.table(fit_id=fit_data$fit_id,outer_fold=fit_data$outer_fold,inner_fold=fit_data$inner_fold,method='elastic_net',endpoint=ends[j],seed=NA_integer_,n_training=length(tr),n_query=length(te),iteration_limit=100000L,reported_convergence=NA_integer_,objective=NA_real_,objective_finite=NA,npasses=as.integer(fit$npasses),jerr=as.integer(fit$jerr),requested_lambda_count=length(lambda_grid),fitted_lambda_count=length(fit$lambda),requested_lambda_min=min(lambda_grid),requested_lambda_max=max(lambda_grid),fitted_lambda_min=min(fit$lambda),fitted_lambda_max=max(fit$lambda),lambda_grid_complete=full_grid,parameters_finite=all(is.finite(fit$beta@x))&&all(is.finite(fit$a0)),status=if(fit$jerr==0L&&full_grid)'completed_requested_lambda_path'else 'solver_or_incomplete_path_warning')

  enet_array[,j,]<-pmin(pmax(as.matrix(predict(fit,xx$test,s=lambda_grid)),0),1)
  if(!is.null(selected_lambda))enet_selected[,j]<-as.numeric(predict(fit,xx$test,s=selected_lambda[[ends[j]]]))
 }
 if(!is.null(selected_lambda))rr$elastic_net<-list(pred=enet_selected,uncertainty=rep(NA_real_,length(te)))
 rows<-rbindlist(lapply(names(rr),function(me){
  pred<-pmin(pmax(rr[[me]]$pred,0),1);unc<-rr[[me]]$uncertainty
  hard<-if(me%in%response_methods)!unseen&nearest<=dist_cut else rep(TRUE,length(te))
  if(me=='kernel_transition')hard<-hard&rr[[me]]$support_n>=10L
  dt<-data.table(fit_id=fit_data$fit_id,outer_fold=fit_data$outer_fold,inner_fold=fit_data$inner_fold,row_id=d$row_id[te],pair_id=d$pair_id[te],cell_line=d$cell_line[te],baseline_context=d$baseline_context[te],perturbation=d$perturbation[te],time_hours=d$time_hours[te],method=me,hard_support=hard,uncertainty=unc,unseen_intervention=unseen,nearest_training_distance=nearest,training_distance_cutoff=dist_cut,n_same_intervention=vapply(d$perturbation[te],function(p)sum(d$perturbation[tr]==p),integer(1)))
  for(j in 1:3){dt[,(paste0('predicted_',ends[j])):=pred[,j]];dt[,(paste0('observed_',ends[j])):=ys[,j]]}
  dt[,mse:=rowMeans((as.matrix(.SD)-ys)^2),.SDcols=paste0('predicted_',ends)]
  dt
 }))
 enet_grid<-rbindlist(lapply(seq_along(lambda_grid),function(l){
  rbindlist(lapply(1:3,function(j)data.table(fit_id=fit_data$fit_id,outer_fold=fit_data$outer_fold,inner_fold=fit_data$inner_fold,row_id=d$row_id[te],pair_id=d$pair_id[te],cell_line=d$cell_line[te],endpoint=ends[j],lambda=lambda_grid[l],prediction=enet_array[,j,l],observed=ys[,j],squared_error=(enet_array[,j,l]-ys[,j])^2)))
 }))
 list(predictions=rows,elastic_grid=enet_grid)
}

all_outer<-list();all_inner<-list();all_lambda<-list();selected_lambda_rows<-list();cutoff_rows<-list();calibrated_rows<-list();extra_outer_rows<-list();extra_cal_rows<-list()
for(of in sort(unique(splits$fold))){
 inner_preds<-list();inner_enet<-list()
 for(inf in 1:4){
  f<-read_latent_fit(paste0('outer_',of,'_inner_',inf));pr<-predict_fit(f)
  inner_preds[[inf]]<-pr$predictions;inner_enet[[inf]]<-pr$elastic_grid
 }
 ip<-rbindlist(inner_preds);eg<-rbindlist(inner_enet)
 stopifnot(!anyDuplicated(ip[,.(row_id,method)]),!anyDuplicated(eg[,.(row_id,endpoint,lambda)]))
 expected_inner<-manifest[outer_fold==of&inner_fold==0L&task_role=='train',row_id]
 stopifnot(setequal(ip$row_id,expected_inner),setequal(eg$row_id,expected_inner))
 scores<-eg[,.(inner_mse=mean(squared_error),n_inner_cases=.N),by=.(outer_fold,endpoint,lambda)]
 selected<-scores[order(inner_mse,-lambda),.SD[1L],by=.(outer_fold,endpoint)]
 selected[,selection_rule:='minimum pooled grouped-inner OOF endpoint MSE; ties prefer stronger penalty']
 ls<-setNames(selected$lambda,selected$endpoint)
 selected_lambda_rows[[length(selected_lambda_rows)+1L]]<-selected
 all_lambda[[length(all_lambda)+1L]]<-eg
 # Freeze all inner-derived cutoffs before computing any outer predictions/errors.
 fold_cuts<-list()
 for(me in response_methods){
  a<-ip[method==me&hard_support&is.finite(uncertainty)]
  taus<-sort(unique(c(0,a$uncertainty,Inf)))
  candidates<-rbindlist(lapply(taus,function(t){b<-a[uncertainty<=t];data.table(threshold=t,n_accepted=nrow(b),n_lines=uniqueN(b$cell_line),inner_rmse=if(nrow(b))sqrt(mean(b$mse))else NA_real_)}))
  for(budget in risk_budgets){
   good<-candidates[n_accepted>=10L&n_lines>=3L&inner_rmse<=budget]
   pick<-if(nrow(good))good[order(-n_accepted,-threshold)][1L]else data.table(threshold=-Inf,n_accepted=0L,n_lines=0L,inner_rmse=NA_real_)
   fold_cuts[[length(fold_cuts)+1L]]<-cbind(data.table(outer_fold=of,method=me,risk_budget=budget,n_inner_total=nrow(ip[method==me]),n_inner_hard_supported=nrow(a)),pick)
  }
 }
 frozen_cuts<-rbindlist(fold_cuts)
 stopifnot(nrow(frozen_cuts)==9L,!anyDuplicated(frozen_cuts[,.(method,risk_budget)]))
 cutoff_rows[[length(cutoff_rows)+1L]]<-copy(frozen_cuts)
 # The outer fit uses the outer-trained scVI space, not any inner coordinates.
 outer_fit<-read_latent_fit(paste0('outer_',of))
 if(include_low_support){
  ex<-read_extra_contexts(paste0('outer_',of),of)
  stopifnot(!length(intersect(outer_fit$data[task_role=='train']$cell_line,ex$cell_line)),all(ex$perturbation%in%perturbation_vocabulary))
  outer_fit$data<-rbindlist(list(outer_fit$data,ex),use.names=TRUE,fill=TRUE)
 }
 op_all<-predict_fit(outer_fit,ls)$predictions
 op<-op_all[row_id<=120L]
 ep<-op_all[row_id>120L]
 op[,accepted:=ifelse(method%in%response_methods,hard_support&is.finite(uncertainty)&uncertainty<=.15,TRUE)]
 op[,selection:='original_heuristic_0.15']
 stopifnot(nrow(op)==24L*6L,all(op$outer_fold==of),all(op$inner_fold==0L))
 all_outer[[length(all_outer)+1L]]<-op
 if(include_low_support){
  stopifnot(nrow(ep)==6L*nrow(extra[outer_fold==of&known_intervention==TRUE]))
  ep[,`:=`(accepted=ifelse(method%in%response_methods,hard_support&is.finite(uncertainty)&uncertainty<=.15,TRUE),selection='original_heuristic_0.15')]
  extra_outer_rows[[length(extra_outer_rows)+1L]]<-ep
 }
 all_inner[[length(all_inner)+1L]]<-ip
 for(k in seq_len(nrow(frozen_cuts))){
  pick<-frozen_cuts[k];me<-pick$method;budget<-pick$risk_budget
  rr<-copy(op[method==me]);rr[,`:=`(risk_budget=budget,selected_uncertainty_threshold=pick$threshold,accepted=hard_support&is.finite(uncertainty)&uncertainty<=pick$threshold,selection='nested_inner_calibrated')]
  calibrated_rows[[length(calibrated_rows)+1L]]<-rr
  if(include_low_support){
   er<-copy(ep[method==me]);er[,`:=`(risk_budget=budget,selected_uncertainty_threshold=pick$threshold,accepted=hard_support&is.finite(uncertainty)&uncertainty<=pick$threshold,selection='nested_inner_calibrated')]
   extra_cal_rows[[length(extra_cal_rows)+1L]]<-er
  }
 }
 cat('Completed outer fold',of,'using separately supplied inner and outer scVI fit tables.\n')
}
outer<-rbindlist(all_outer);inner<-rbindlist(all_inner);grid<-rbindlist(all_lambda);selected<-rbindlist(selected_lambda_rows);cuts<-rbindlist(cutoff_rows);cal<-rbindlist(calibrated_rows)
stopifnot(nrow(outer)==720L,!anyDuplicated(outer[,.(row_id,method)]),nrow(cal)==1080L,!anyDuplicated(cal[,.(row_id,method,risk_budget)]),nrow(selected)==15L,nrow(cuts)==45L)
optimizer_audit<-rbindlist(optimizer_rows)
stopifnot(nrow(optimizer_audit)==150L,uniqueN(optimizer_audit$fit_id)==25L,!anyDuplicated(optimizer_audit[,.(fit_id,method,endpoint,seed)]))
write_tsv(optimizer_audit,'nested_scvi_optimizer_audit.tsv')
write_tsv(optimizer_audit[,.(n_fits=.N,all_parameters_finite=all(parameters_finite)),by=.(method,status)],'nested_scvi_optimizer_status_counts.tsv')
write_tsv(rbindlist(fit_table_checks),'nested_scvi_fit_table_checks.tsv')
write_tsv(outer,'nested_scvi_outer_predictions.tsv');write_tsv(inner,'nested_scvi_inner_response_predictions.tsv');write_tsv(grid,'nested_scvi_inner_elastic_lambda_predictions.tsv');write_tsv(selected,'nested_scvi_selected_lambdas.tsv');write_tsv(cuts,'nested_scvi_selected_thresholds.tsv');write_tsv(cal,'nested_scvi_calibrated_outer_predictions.tsv')

# Conditional fixed-decision uncertainty: no scVI/model/cutoff refitting here.
B<-2000L;set.seed(2026110308L)
lines<-sort(unique(tasks$cell_line));line_draws<-rmultinom(B,size=length(lines),prob=rep(1/length(lines),length(lines)))
percentile<-function(v){v<-v[is.finite(v)];if(length(v))unname(quantile(v,c(.025,.975)))else c(NA_real_,NA_real_)}
summarize_model<-function(a,selection,risk_budget=NA_real_,elastic_reference=outer,expected_n=120L,cluster_lines=lines,cluster_draws=line_draws){
 a<-copy(a);stopifnot(nrow(a)==expected_n,!anyDuplicated(a$row_id));el<-elastic_reference[method=='elastic_net'][match(a$row_id,row_id)]
 stopifnot(identical(a$row_id,el$row_id),identical(a$pair_id,el$pair_id))
 w<-cluster_draws[match(a$cell_line,cluster_lines),,drop=FALSE];mask<-as.numeric(a$accepted);den<-colSums(w*mask);tot<-colSums(w);coverage_ci<-percentile(den/tot);n_lines<-uniqueN(a[accepted==TRUE]$cell_line)
 rbindlist(lapply(c('overall',ends),function(endpoint){
  err<-if(endpoint=='overall')a$mse else (a[[paste0('predicted_',endpoint)]]-a[[paste0('observed_',endpoint)]])^2
  eerr<-if(endpoint=='overall')el$mse else (el[[paste0('predicted_',endpoint)]]-el[[paste0('observed_',endpoint)]])^2
  stopifnot(all(is.finite(err[a$accepted])),all(is.finite(eerr[a$accepted])))
  n<-sum(mask);rmse<-if(n)sqrt(mean(err[a$accepted]))else NA_real_;ermse<-if(n)sqrt(mean(eerr[a$accepted]))else NA_real_
  masked_err<-ifelse(a$accepted,err,0);masked_elastic_err<-ifelse(a$accepted,eerr,0)
  rb<-if(n_lines>=2L)sqrt(colSums(w*masked_err)/den)else rep(NA_real_,B)
  eb<-if(n_lines>=2L)sqrt(colSums(w*masked_elastic_err)/den)else rep(NA_real_,B)
  ci<-percentile(rb);eci<-percentile(eb);dci<-percentile(rb-eb)
  data.table(selection=selection,risk_budget=risk_budget,method=unique(a$method),endpoint=endpoint,n_total=nrow(a),n_accepted=n,n_accepted_cell_lines=n_lines,coverage=n/nrow(a),coverage_ci_low=coverage_ci[1],coverage_ci_high=coverage_ci[2],rmse=rmse,rmse_ci_low=ci[1],rmse_ci_high=ci[2],matched_elastic_rmse=ermse,matched_elastic_ci_low=eci[1],matched_elastic_ci_high=eci[2],difference_model_minus_elastic=rmse-ermse,difference_ci_low=dci[1],difference_ci_high=dci[2],valid_rmse_bootstraps=sum(is.finite(rb)),uncertainty_scope='cell-line cluster bootstrap of fixed OOF decisions/predictions; no pipeline refit')
 }))
}
heuristic_summary<-rbindlist(lapply(all_methods,function(me)summarize_model(outer[method==me],'original_heuristic_0.15')))
cal_summary<-rbindlist(lapply(response_methods,function(me)rbindlist(lapply(risk_budgets,function(b)summarize_model(cal[method==me&risk_budget==b],'nested_inner_calibrated',b)))))
write_tsv(heuristic_summary,'nested_scvi_heuristic_metrics.tsv');write_tsv(cal_summary,'nested_scvi_calibrated_metrics.tsv')
write_tsv(outer[,.(n_total=.N,n_accepted=sum(accepted),rmse=if(any(accepted))sqrt(mean(mse[accepted]))else NA_real_),by=.(outer_fold,method)],'nested_scvi_heuristic_by_fold.tsv')
write_tsv(cal[,.(n_total=.N,n_accepted=sum(accepted),rmse=if(any(accepted))sqrt(mean(mse[accepted]))else NA_real_),by=.(outer_fold,method,risk_budget)],'nested_scvi_calibrated_by_fold.tsv')
if(include_low_support){
 extra_outer<-rbindlist(extra_outer_rows);extra_cal<-rbindlist(extra_cal_rows)
 stopifnot(nrow(extra_outer)==144L*6L,!anyDuplicated(extra_outer[,.(row_id,method)]),nrow(extra_cal)==144L*9L,!anyDuplicated(extra_cal[,.(row_id,method,risk_budget)]))
 extra_scope[,context_present:=TRUE] # All 210 context joins asserted before excluding novel identities.
 enrich_extra<-function(a)merge(a,extra_scope[,.(row_id,baseline_cells,treated_cells,minimum_support,maximum_program_se,count_bin,se_bin,support_deficit)],by='row_id',all.x=TRUE,sort=FALSE)
 extra_outer<-enrich_extra(extra_outer);extra_cal<-enrich_extra(extra_cal)
 strata<-rbindlist(list(data.table(stratum_type='all_known_identities',stratum='all',row_id=extra[known_intervention==TRUE,row_id]),extra[known_intervention==TRUE,.(stratum_type='minimum_arm_cell_count',stratum=count_bin,row_id)],extra[known_intervention==TRUE,.(stratum_type='maximum_program_SE',stratum=se_bin,row_id)]))
 extra_strata<-unique(strata[,.(stratum_type,stratum)])
 set.seed(2026110309L)
 extra_bootstraps<-lapply(seq_len(nrow(extra_strata)),function(i){ss<-extra_strata[i];ids<-strata[stratum_type==ss$stratum_type&stratum==ss$stratum,row_id];ll<-sort(unique(extra[row_id%in%ids]$cell_line));list(lines=ll,draws=rmultinom(B,size=length(ll),prob=rep(1/length(ll),length(ll))))})
 summarize_extra<-function(a,sel,budget=NA_real_){
  rbindlist(lapply(seq_len(nrow(extra_strata)),function(i){
   ss<-extra_strata[i];ids<-strata[stratum_type==ss$stratum_type&stratum==ss$stratum,row_id]
   bb<-extra_bootstraps[[i]];v<-summarize_model(a[row_id%in%ids],sel,budget,extra_outer,length(ids),bb$lines,bb$draws)
   cbind(ss,data.table(evaluation_panel='additional low-support conditions; same benchmark lines and known intervention vocabulary',n_total_cell_lines=uniqueN(a[row_id%in%ids]$cell_line)),v)
  }))
 }
 extra_hm<-rbindlist(lapply(all_methods,function(me)summarize_extra(extra_outer[method==me],'original_heuristic_0.15')))
 extra_cm<-rbindlist(lapply(response_methods,function(me)rbindlist(lapply(risk_budgets,function(b)summarize_extra(extra_cal[method==me&risk_budget==b],'nested_inner_calibrated',b)))))
 write_tsv(extra_scope,'nested_scvi_low_support_scope.tsv')
 write_tsv(extra_scope[,.(n_conditions=.N,n_cell_lines=uniqueN(cell_line)),by=.(known_intervention,perturbation,support_deficit)],'nested_scvi_low_support_identity_counts.tsv')
 write_tsv(extra_outer,'nested_scvi_low_support_predictions.tsv');write_tsv(extra_cal,'nested_scvi_low_support_calibrated_predictions.tsv')
 write_tsv(extra_hm,'nested_scvi_low_support_heuristic_metrics.tsv');write_tsv(extra_cm,'nested_scvi_low_support_calibrated_metrics.tsv')
 write_tsv(extra_outer[,.(n_total=.N,n_accepted=sum(accepted),rmse=if(any(accepted))sqrt(mean(mse[accepted]))else NA_real_),by=.(outer_fold,method)],'nested_scvi_low_support_heuristic_by_fold.tsv')
 write_tsv(extra_cal[,.(n_total=.N,n_accepted=sum(accepted),rmse=if(any(accepted))sqrt(mean(mse[accepted]))else NA_real_),by=.(outer_fold,method,risk_budget)],'nested_scvi_low_support_calibrated_by_fold.tsv')
}
checks<-data.table(check=c('25 distinct scVI fit identifiers','five outer folds','120 test tasks each with six models','1080 calibrated decisions','15 inner-selected endpoint lambdas','45 training-only cutoff choices','zero-acceptance metrics represented','150 optimizer statuses recorded without extra fits','all 75 EN fits supply the complete fixed 17-penalty path'),passed=c(uniqueN(manifest$fit_id)==25L,uniqueN(outer$outer_fold)==5L,nrow(outer)==720L&&!anyDuplicated(outer[,.(row_id,method)]),nrow(cal)==1080L&&!anyDuplicated(cal[,.(row_id,method,risk_budget)]),nrow(selected)==15L,nrow(cuts)==45L,all(is.na(cal_summary[n_accepted==0]$rmse)),nrow(optimizer_audit)==150L&&all(optimizer_audit$parameters_finite),all(optimizer_audit[method=='elastic_net']$jerr==0L&optimizer_audit[method=='elastic_net']$lambda_grid_complete&optimizer_audit[method=='elastic_net']$fitted_lambda_count==17L)))
if(include_low_support)checks<-rbind(checks,data.table(check=c('210 extra conditions have outer-held-line contexts','144 known-vocabulary extra conditions each with six models','66 novel-identity conditions separately audited without vocabulary expansion','1296 separately calibrated extra decisions'),passed=c(nrow(extra_scope)==210L&&all(extra_scope$context_present),nrow(extra_outer)==864L&&!anyDuplicated(extra_outer[,.(row_id,method)]),sum(!extra_scope$known_intervention)==66L&&all(extra_scope[known_intervention==FALSE]$prespecified_response_identity_rejection),nrow(extra_cal)==1296L&&!anyDuplicated(extra_cal[,.(row_id,method,risk_budget)]))))
write_tsv(checks,'nested_scvi_downstream_verification.tsv');stopifnot(all(checks$passed))
writeLines(capture.output(sessionInfo()),file.path(out,'nested_scvi_R_sessionInfo.txt'))
cat('Strict nested scVI downstream model selection and evaluation completed.\n')

invisible(list(outer_predictions=outer,calibrated_predictions=cal,heuristic_metrics=heuristic_summary,calibrated_metrics=cal_summary,selected_lambdas=selected,selected_thresholds=cuts,low_support_predictions=if(include_low_support)extra_outer else NULL,low_support_calibrated_predictions=if(include_low_support)extra_cal else NULL))
}

# Compare two prespecified budgets and export the NN sensitivity without duplicate baselines.
compare_nested_scvi_nn_budgets<-function(primary_dir,sensitivity_dir,out,primary_result=NULL,sensitivity_result=NULL){
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 read_result<-function(d,n)fread(file.path(d,paste0('nested_scvi_',n,'.tsv')))
 nn<-'conditional_rate_network'
 checks<-list()
 same_vectors<-function(a,b,label){
  keys<-intersect(c('fit_id','outer_fold','inner_fold','row_id','method','risk_budget','endpoint','lambda'),names(a))
  if(length(keys)){setorderv(a,keys);setorderv(b,keys)}
  ok<-identical(names(a),names(b))&&nrow(a)==nrow(b)&&all(vapply(names(a),function(z)identical(a[[z]],b[[z]]),logical(1)))
  checks[[length(checks)+1L]]<<-data.table(check=label,passed=ok)
  if(!ok)stop('Budget invariance failed: ',label)
 }
 if(!is.null(primary_result)&&!is.null(sensitivity_result)){
  for(n in c('outer_predictions','calibrated_predictions','low_support_predictions','low_support_calibrated_predictions'))same_vectors(copy(primary_result[[n]][method!=nn]),copy(sensitivity_result[[n]][method!=nn]),paste0('in-memory non-NN ',n,' bitwise identical'))
  same_vectors(copy(primary_result$selected_lambdas),copy(sensitivity_result$selected_lambdas),'in-memory EN lambdas bitwise identical')
  same_vectors(copy(primary_result$selected_thresholds[method!=nn]),copy(sensitivity_result$selected_thresholds[method!=nn]),'in-memory non-NN cutoffs bitwise identical')
 }
 for(n in c('outer_predictions','calibrated_outer_predictions','inner_response_predictions','low_support_predictions','low_support_calibrated_predictions')){
  a<-read_result(primary_dir,n);b<-read_result(sensitivity_dir,n)
  same_vectors(a[method!=nn],b[method!=nn],paste0('non-NN ',n,' identical'))
 }
 same_vectors(read_result(primary_dir,'selected_lambdas'),read_result(sensitivity_dir,'selected_lambdas'),'all 15 EN lambdas identical')
 same_vectors(read_result(primary_dir,'inner_elastic_lambda_predictions'),read_result(sensitivity_dir,'inner_elastic_lambda_predictions'),'all inner EN grid predictions identical')
 a<-read_result(primary_dir,'selected_thresholds');b<-read_result(sensitivity_dir,'selected_thresholds')
 same_vectors(a[method!=nn],b[method!=nn],'all 30 non-NN inner cutoffs identical')
 aa<-read_result(primary_dir,'optimizer_audit')[method==nn];bb<-read_result(sensitivity_dir,'optimizer_audit')[method==nn]
 stopifnot(nrow(aa)==75L,nrow(bb)==75L,all(aa$iteration_limit==600L),all(bb$iteration_limit==6000L))
 audit<-merge(aa,bb,by=c('fit_id','outer_fold','inner_fold','method','seed'),suffixes=c('_600','_6000'))
 stopifnot(nrow(audit)==75L,all(audit$n_training_600==audit$n_training_6000),all(audit$n_query_600==audit$n_query_6000))
 audit[,`:=`(objective_reduction=objective_600-objective_6000,relative_objective_reduction=(objective_600-objective_6000)/pmax(abs(objective_600),.Machine$double.eps))]
 fwrite(audit,file.path(out,'nested_scvi_nn_optimizer_budget_comparison.tsv'),sep='\t')
 fwrite(bb,file.path(out,'nested_scvi_nn_6000_optimizer_audit.tsv'),sep='\t')
 for(n in c('outer_predictions','calibrated_outer_predictions','inner_response_predictions','selected_thresholds','heuristic_metrics','calibrated_metrics','heuristic_by_fold','calibrated_by_fold','low_support_predictions','low_support_calibrated_predictions','low_support_heuristic_metrics','low_support_calibrated_metrics','low_support_heuristic_by_fold','low_support_calibrated_by_fold')){
  a<-read_result(sensitivity_dir,n)[method==nn]
  a[,nn_iteration_limit:=6000L]
  fwrite(a,file.path(out,paste0('nested_scvi_nn_6000_',n,'.tsv')),sep='\t')
 }
 cv<-rbindlist(checks);fwrite(cv,file.path(out,'nested_scvi_nn_budget_invariance.tsv'),sep='\t')
 invisible(list(checks=cv,optimizer_comparison=audit))
}

# Five compact direct-inspection tables; original published-result tables are untouched.
export_nested_scvi_tables<-function(bundle,primary_dir,sensitivity_dir,compact_dir,out){
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 inp<-readRDS(file.path(bundle,'inputs/mixseq_nested_scvi.rds'))
 rr<-function(d,n)fread(file.path(d,paste0('nested_scvi_',n,'.tsv')))
 p<-list();m<-list()
 for(q in c('primary','additional_known_vocabulary')){
  prefix<-if(q=='primary')''else 'low_support_'
  for(policy in c('heuristic','calibrated')){
   pred_name<-if(policy=='heuristic')paste0(prefix,if(q=='primary')'outer_predictions'else 'predictions')else paste0(prefix,if(q=='primary')'calibrated_outer_predictions'else 'calibrated_predictions')
   # Additional rows use low_support_predictions / low_support_calibrated_predictions.
   for(budget in c(600L,6000L)){
    d<-if(budget==600L)primary_dir else sensitivity_dir
    a<-rr(d,pred_name);v<-rr(d,paste0(prefix,policy,'_metrics'))
    if(budget==6000L){a<-a[method=='conditional_rate_network'];v<-v[method=='conditional_rate_network']}
    a[,`:=`(evaluation_panel=q,nn_iteration_limit=budget)];v[,`:=`(evaluation_panel=q,nn_iteration_limit=budget)]
    if(!'stratum_type'%in%names(v))v[,`:=`(stratum_type='all_primary',stratum='all')]
    p[[length(p)+1L]]<-a;m[[length(m)+1L]]<-v
   }
  }
 }
 predictions<-rbindlist(p,use.names=TRUE,fill=TRUE);metrics<-rbindlist(m,use.names=TRUE,fill=TRUE)
 all_tasks<-as.data.table(readRDS(file.path(bundle,'inputs/mixseq_modeling.rds')))
 counts<-all_tasks[,.(pair_id,baseline_cells,treated_cells,minimum_support=pmin(baseline_cells,treated_cells),maximum_program_se=pmax(p_apoptosis_se,p_proliferation_se,p_stress_se))]
 predictions[,c(intersect(setdiff(names(counts),'pair_id'),names(predictions))):=NULL]
 predictions<-merge(predictions,counts,by='pair_id',sort=FALSE,all.x=TRUE)
 scope<-rr(primary_dir,'low_support_scope')
 man<-rr(primary_dir,'fold_manifest')[inner_fold==0L]
 main<-unique(man[,.(row_id,pair_id,cell_line,baseline_context,perturbation)])
 main<-merge(main,counts,by='pair_id',sort=FALSE,all.x=TRUE)
 sf<-man[task_role=='test',.(pair_id,outer_fold)]
 main<-merge(main,sf,by='pair_id',sort=FALSE)
 main[,`:=`(known_intervention=TRUE,numerically_evaluated=TRUE,prespecified_response_identity_rejection=FALSE,context_present=TRUE,evaluation_scope='original adequate-support 120-condition primary panel',evaluation_panel='primary')]
 scope[,evaluation_panel:=ifelse(known_intervention,'additional_known_vocabulary','additional_novel_identity')]
 eligibility<-rbindlist(list(main,scope),use.names=TRUE,fill=TRUE)
 train<-copy(as.data.table(inp$fit_audit));train[,`:=`(record_type='representation_scVI',method='scVI',nn_iteration_limit=NA_integer_)]
 op<-rr(primary_dir,'optimizer_audit');op[,`:=`(record_type='prediction_optimizer',nn_iteration_limit=600L)]
 os<-rr(sensitivity_dir,'optimizer_audit')[method=='conditional_rate_network'];os[,`:=`(record_type='prediction_optimizer',nn_iteration_limit=6000L)]
 training<-rbindlist(list(train,op,os),use.names=TRUE,fill=TRUE)
 lam<-rr(primary_dir,'selected_lambdas');lam[,`:=`(record_type='endpoint_lambda',method='elastic_net',nn_iteration_limit=600L)]
 cu<-rr(primary_dir,'selected_thresholds');cu[,`:=`(record_type='inner_uncertainty_cutoff',nn_iteration_limit=600L)]
 cs<-rr(sensitivity_dir,'selected_thresholds')[method=='conditional_rate_network'];cs[,`:=`(record_type='inner_uncertainty_cutoff',nn_iteration_limit=6000L)]
 policies<-rbindlist(list(lam,cu,cs),use.names=TRUE,fill=TRUE)
 tables<-list(predictions=predictions,metrics=metrics,eligibility=eligibility,training=training,policies=policies)
 for(n in names(tables)){
  fwrite(tables[[n]],file.path(out,paste0(n,'.tsv')),sep='\t')
  con<-gzfile(file.path(bundle,'tables',paste0('mixseq_nested_',n,'.tsv.gz')),'wt',compression=9)
  write.table(tables[[n]],con,sep='\t',quote=TRUE,row.names=FALSE,na='NA');close(con)
 }
 stopifnot(nrow(predictions)==5016L,nrow(metrics)==456L,nrow(eligibility)==330L,nrow(training)==250L,nrow(policies)==75L)
 invisible(tables)
}

if(sys.nframe()==0L){
 args<-commandArgs(trailingOnly=TRUE)
 if(length(args)<3L)stop("Usage: Rscript nested_scvi_models.R BUNDLE_ROOT LATENT_ROOT_OR_RDS OUTPUT_DIR")
 extra_args<-if(length(args)>3L)args[-(1:3)]else character()
 nn_arg<-grep('^--nn-maxit=',extra_args,value=TRUE)
 if(length(nn_arg)>1L||any(!extra_args%in%c('--primary-only',nn_arg)))stop('Unknown or duplicate optional arguments.')
 nn_limit<-if(length(nn_arg))as.integer(sub('^--nn-maxit=','',nn_arg))else 600L
 run_nested_scvi_models(args[1],args[2],args[3],include_low_support=!'--primary-only'%in%extra_args,nn_maxit=nn_limit)
}
