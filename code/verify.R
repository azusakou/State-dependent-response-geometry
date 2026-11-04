suppressPackageStartupMessages(library(data.table))
r<-readRDS('reference/published_numerical_reference.rds');checks<-list()
check<-function(name,value,limit=1e-8){ok<-if(is.logical(value))all(value)else is.finite(value)&&value<=limit;checks[[length(checks)+1]]<<-data.table(check=name,passed=ok,value=as.character(value),tolerance=as.character(limit));if(!ok)warning('Numerical check failed: ',name)}
new<-fread('results/task_registry/state_transition_tasks.tsv');old<-copy(r$tasks);setorder(new,task_id);setorder(old,task_id)
check('complete task identifiers and eligibility',identical(new$task_id,old$task_id)&&identical(new$status,old$status));check('715 eligible tasks',sum(new$status=='eligible')==715)
check('Norman basin assignment',all(new$target_state_id==old$target_state_id|is.na(new$target_state_id)&is.na(old$target_state_id)))
n<-fread('results/norman_control/combination_interactions.tsv');c<-merge(n,r$norman_pairs,by='combination');check('Norman primary projection',max(abs(c$interaction_projection.x-c$interaction_projection.y)));check('Norman 131 classes',all(c$interaction_class.x==c$interaction_class.y));check('Norman technical-GEM p',max(abs(c$replicate_projection_p.x-c$replicate_projection_p.y)))
m<-fread('results/mixseq_control_models/metrics.tsv');z<-merge(m,r$mixseq_metrics,by=c('split','method'));check('MIX model coverage',max(abs(z$coverage.x-z$coverage.y)));check('MIX model RMSE',max(abs(z$target_state_rmse.x-z$target_state_rmse.y),na.rm=TRUE),1e-6)
s<-fread('results/mixseq_control_models/split_assignments.tsv');o<-as.data.table(r$mixseq_splits);setorder(s,split,row_id);setorder(o,split,row_id);check('MIX atomic split membership',identical(s,o))
for(path in c('norman_metrics','norman_baseline_metrics')){s<-fread(if(path=='norman_metrics')'results/grouped_validation/norman_metrics.tsv'else'results/extended_baselines/norman_baseline_metrics.tsv');z<-merge(s,r[[path]],by=c('split','method'));check(paste(path,'interaction RMSE'),max(abs(z$interaction_rmse.x-z$interaction_rmse.y)),1e-5);check(paste(path,'AUROC'),max(abs(z$synergy_auc.x-z$synergy_auc.y),na.rm=TRUE),1e-6)}
e<-fread('results/control_atlas/control_node_evidence.tsv');z<-merge(e,r$evidence,by='gene');check('9894 gene evidence classes',nrow(z)==9894&&all(z$evidence_class.x==z$evidence_class.y))
a<-fread('results/multiperturb_multiome/target_reproducibility.tsv');z<-merge(a,r$multiome,by=c('context','guide_gene'));check('multiome joint guide criteria',all(z$passes_task_qc.x==z$passes_task_qc.y));check('multiome RNA guide cosine',max(abs(z$median_rna_guide_cosine.x-z$median_rna_guide_cosine.y),na.rm=TRUE))
s<-fread('results/control_atlas/observed_support_sequences.tsv');z<-merge(s,r$sequence,by=c('cell_line','perturbation'));check('retrospective exposure labels',nrow(z)==158&&all(z$observed_support_policy.x==z$observed_support_policy.y))
s<-fread('results/control_atlas/control_degeneracy_sensitivity.tsv');z<-merge(s,r$degeneracy,by=c('target_state_id','cosine_threshold'));check('endpoint diversity across 3 thresholds',max(abs(z$control_degeneracy.x-z$control_degeneracy.y)))
if(file.exists('results/lineage_sensitivity/lineage_sensitivity.tsv')){
 s<-fread('results/lineage_sensitivity/lineage_sensitivity.tsv');o<-as.data.table(r$lineage_sensitivity);z<-merge(s,o,by=c('dataset','setting','adjustment'));check('LARRY stress sensitivity coefficients',max(abs(z$beta.x-z$beta.y)))
 s<-fread('results/lineage_sensitivity/lineage_full_pipeline_bootstrap_summary.tsv');z<-merge(s,r$lineage_bootstrap,by='dataset');check('LARRY full-pipeline bootstrap interval',max(abs(z$lower.x-z$lower.y),abs(z$upper.x-z$upper.y)))
}
# Independent-replay references for the added review checks.
r2<-readRDS('reference/additional_numerical_reference.rds')
for(path in names(r2)) {
 s<-fread(file.path('results',path));o<-as.data.table(r2[[path]])
 ok<-identical(names(s),names(o))&&nrow(s)==nrow(o);delta<-0
 if(ok) for(nm in names(s)) {
  a<-s[[nm]];b<-o[[nm]]
  if(!identical(is.na(a),is.na(b))){ok<-FALSE;break}
  if(is.numeric(a)&&is.numeric(b)) {
   keep<-is.finite(a)&is.finite(b)
   if(any(keep))delta<-max(delta,max(abs(a[keep]-b[keep])))
   if(!identical(a[!is.finite(a)],b[!is.finite(b)]))ok<-FALSE
  } else if(!identical(as.character(a),as.character(b)))ok<-FALSE
 }
 check(paste('numerical replay',path),if(ok)delta else Inf,1e-6)
}
n<-fread('results/prediction_validation/nested_calibration_outer_predictions.tsv')
check('nested outer decisions complete and unique',nrow(n)==1080L&&!anyDuplicated(n[,.(row_id,method,risk_budget)]))
m<-fread('results/prediction_validation/nested_calibration_membership.tsv')
check('nested calibration excludes outer cell lines',all(m[,.(ok=length(intersect(cell_line[role=='outer_train'],cell_line[role=='outer_test']))==0),by=outer_fold]$ok))
m<-fread('results/prediction_validation/rawRNA_PCA_membership.tsv')
check('RNA training representation excludes test contexts',all(m[,.(ok=length(intersect(baseline_context[role=='fit_features_PCA'],baseline_context[role=='test_project_only']))==0),by=outer_fold]$ok))
bm<-fread('tables/baseline_cell_context_membership.tsv.gz');bp<-readRDS('inputs/mixseq_raw_baseline_profiles.rds')
check('7328 baseline cells map to 57 source contexts',nrow(bm)==7328L&&uniqueN(bm$baseline_context)==57L&&!anyDuplicated(bm[,.(source_file,cell_id)])&&sum(bp$metadata$n_baseline)==nrow(bm))

g<-fread('results/prediction_validation/functional_match_summary.tsv')
check('orthogonal growth matching and coverage',g$all_pairs==131L&&g$matched_pairs==121L&&g$eligible_pairs==117L&&g$matched_eligible_pairs==110L&&g$matched_eligible_basins==15L)
# Fully nested learned-representation replay and honest fixed-budget diagnostics.
ns<-readRDS('inputs/mixseq_nested_scvi.rds');nv<-ns$verification
check('nested scVI completed25 slots and15 distinct training states',isTRUE(nv$complete)&&nv$fits==25L&&nv$distinct_final_encoder_states==15L&&nv$distinct_training_cell_sets==15L&&isTRUE(nv$identical_state_only_with_identical_training_membership))
check('nested scVI fresh fits and frozen query encoding',isTRUE(nv$all_scvi_fits_fresh)&&isTRUE(nv$all_frozen_query_weight_hashes_identical)&&isTRUE(nv$all_training_batches_present))
fa<-as.data.table(ns$fit_audit);ge<-as.data.table(ns$selected_genes);eh<-as.data.table(ns$epoch_metrics);lz<-as.data.table(ns$latents)
check('nested scVI complete fitted-input audit',nrow(fa)==25L&&nrow(ge)==100000L&&nrow(eh)==1025L&&nrow(lz)==2520L&&nrow(ns$context_latents)==1675L)
check('nested scVI reported selected genes and200epoch budget',all(ge[,.(ok=.N==4000L&&!anyDuplicated(gene)),by=fit_id]$ok)&&all(fa$selected_genes==4000L)&all(fa$actual_epochs==200L)&uniqueN(fa$weight_hash_after_query)==15L&uniqueN(fa$training_cell_order_sha256)==15L)
check('nested scVI finite reported training and scheduled validation histories',all(is.finite(eh$train_loss))&&all(is.finite(eh$elbo_train))&&all(is.finite(eh[epoch%%5L==0L]$validation_loss))&&all(is.finite(eh[epoch%%5L==0L]$elbo_validation))&&all(eh[,.(ok=max(epoch)==200L&&.N==41L),by=fit_id]$ok))
check('nested scVI query leaves all encoder weights fixed',all(fa$weight_hash_before_query==fa$weight_hash_after_query)&&all(tolower(as.character(fa$query_updates_weights))=='false')&&all(fa$query_unseen_batches==0L))
fm<-fread('results/nested_scvi_600/nested_scvi_fold_manifest.tsv');expect<-fm[task_role!='outer_test_excluded'];zz<-merge(lz,expect,by=c('fit_id','row_id','pair_id'),suffixes=c('_input','_expected'))
check('nested scVI exact2520 task memberships and roles',nrow(zz)==2520L&&!anyDuplicated(zz[,.(fit_id,row_id)])&&all(zz$task_role_input==zz$task_role_expected)&&all(zz$cell_line_input==zz$cell_line_expected)&&all(zz$baseline_context_input==zz$baseline_context_expected))
cm<-fread('inputs/mixseq_baseline_cells.tsv.gz');bench<-unique(fm$cell_line);aux<-setdiff(unique(cm$cell_line),bench)
cell_checks<-lapply(seq_len(nrow(fa)),function(i){a<-fa[i];key<-a$fit_id;f<-fm[fit_id==key];ol<-unique(fm[outer_fold==a$outer_fold&inner_fold==0L&task_role=='test']$cell_line);ql<-unique(f[task_role%in%c('test','validation')]$cell_line);tl<-union(aux,unique(f[task_role=='train']$cell_line));length(intersect(tl,ql))==0L&&length(intersect(tl,ol))==0L&&sum(cm$cell_line%in%tl)==a$training_cells&&sum(cm$cell_line%in%ql)==a$query_cells&&length(tl)==a$reference_lines&&length(ql)==a$query_lines&&(if(a$inner_fold==0L)0L else sum(cm$cell_line%in%ol))==a$excluded_outer_cells})
check('nested scVI all25 actual cell memberships agree with frozen line exclusions',nrow(cm)==21936L&&uniqueN(cm$cell_line)==169L&&length(aux)==128L&&all(unlist(cell_checks)))
check('nested scVI no inner validation or outer lines enter model training',all(fm[task_role!='outer_test_excluded',.(ok=length(intersect(cell_line[task_role=='train'],cell_line[task_role!='train']))==0L),by=fit_id]$ok))
pr<-fread('results/nested_scvi/predictions.tsv');ee<-fread('results/nested_scvi/eligibility.tsv');po<-fread('results/nested_scvi/policies.tsv');tt<-fread('results/nested_scvi/training.tsv')
check('nested scVI primary plus additional predictions complete without duplicates',nrow(pr)==5016L&&!anyDuplicated(pr[,.(evaluation_panel,row_id,method,selection,risk_budget,nn_iteration_limit)]))
check('nested scVI eligibility330 rows keeps120 primary144known66novel separate',nrow(ee)==330L&&sum(ee$evaluation_panel=='primary')==120L&&sum(ee$evaluation_panel=='additional_known_vocabulary')==144L&&sum(ee$evaluation_panel=='additional_novel_identity')==66L&&all(ee$context_present)&&!any(ee[evaluation_panel=='additional_novel_identity']$numerically_evaluated))
check('nested scVI policy and audit counts',nrow(po)==75L&&nrow(tt)==250L)
for(budget in c(600L,6000L)){
 oa<-fread(file.path(paste0('results/nested_scvi_',budget),'nested_scvi_optimizer_audit.tsv'))
 en<-oa[method=='elastic_net'];nn<-oa[method=='conditional_rate_network']
 check(paste('nested EN complete17lambda path budget',budget),nrow(en)==75L&&all(en$jerr==0L&en$lambda_grid_complete&en$fitted_lambda_count==17L&en$fitted_lambda_min==1e-4&en$fitted_lambda_max==1))
 check(paste('nested NN finite objectives and recorded status budget',budget),nrow(nn)==75L&&all(nn$objective_finite&nn$parameters_finite)&&all(nn$reported_convergence==if(budget==600L)1L else 0L))
}
bi<-fread('results/nested_scvi_nn6000/nested_scvi_nn_budget_invariance.tsv');check('nested NN budget leaves all nonNN predictions and policies invariant',nrow(bi)==14L&&all(bi$passed))

ch<-rbindlist(checks);fwrite(ch,'results/verification.tsv',sep='\t');print(ch);if(!all(ch$passed))stop('One or more numerical reproduction checks failed; see results/verification.tsv.')
