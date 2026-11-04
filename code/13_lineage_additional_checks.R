#!/usr/bin/env Rscript
# Additional sensitivity sampling and measured-quality diagnostics.
# All estimates are observational; overlap weighting is a descriptive comparison.
suppressPackageStartupMessages({library(data.table);library(BiocNeighbors)})
args<-commandArgs(trailingOnly=TRUE)
pkg<-if(length(args)&&dir.exists(file.path(args[1],'inputs')))args[1]else'.'
out<-if(length(args)>1&&!grepl('^--',args[2]))args[2]else'results/lineage_validation'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache<-readRDS(file.path(pkg,'inputs','lineage_inputs.rds'))
all_inputs<-cache$data;desired_features<-cache$features
base_cov<-c('cell_cycle_state','stress_state','chromatin_regulator_proxy','metabolic_state','mitochondrial_fraction')
stress_axes<-c('hallmark__HYPOXIA','hallmark__REACTIVE_OXYGEN_SPECIES_PATHWAY','hallmark__UNFOLDED_PROTEIN_RESPONSE')
for(e in parse(file.path(pkg,'code','10_lineage.R'))) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('fit_one','neighbours','geometry'))eval(e)
set.seed(202611033)
weighted_average<-function(x,w)sum(x*w)/sum(w)
weighted_var<-function(x,w){mu<-weighted_average(x,w);sum(w*(x-mu)^2)/sum(w)}
cluster_ci<-function(y,high,w,cluster){
 dd<-data.frame(y=y,high=high,w=w,cluster=cluster);f<-lm(y~high,data=dd,weights=w)
 x<-model.matrix(f);e<-resid(f);bread<-solve(crossprod(x*w,x));score<-rowsum(x*(w*e),cluster,reorder=FALSE);g<-nrow(score);n<-length(y);p<-ncol(x)
 vc<-(g/(g-1))*((n-1)/(n-p))*bread%*%crossprod(score)%*%bread
 est<-coef(f)[['high']];se<-sqrt(vc[2,2]);c(estimate=est,lower=est-qt(.975,g-1)*se,upper=est+qt(.975,g-1)*se)
}
balances<-list();overlaps<-list();ps_rows<-list();selection<-list();subset_slopes<-list();condition_slopes<-list();donor_results<-list()
for(nm in names(all_inputs)){
 d<-copy(all_inputs[[nm]]);d[,row_id:=.I]
 # One clone gets total unit weight one even if sampled under several conditions.
 d[,clone_weight:=1/.N,by=clone_id]
 threshold<-median(unique(d,by='clone_id')$stress_state)
 d[,high:=as.integer(stress_state>threshold)]
 y<-d$local_effective_rank_stress_state_loo
 for(spec in c('quality_only','quality_and_other_states')){
   rhs<-c('factor(baseline_n)','splines::ns(log_descendant_n,df=2)','splines::ns(log_detected_genes,df=3)','splines::ns(mitochondrial_fraction,df=3)')
   if(spec=='quality_and_other_states')rhs<-c(rhs,'cell_cycle_state','chromatin_regulator_proxy','metabolic_state')
   if(uniqueN(d$future_condition)>1)rhs<-c(rhs,'factor(future_condition)')
   if(uniqueN(d$dominant_library)>1)rhs<-c(rhs,'factor(dominant_library)')
   # quasibinomial avoids a fractional-prior-weight warning; point estimates are
   # the same weighted-logistic score solution used for overlap balance.
   f<-glm(as.formula(paste('high~',paste(rhs,collapse='+'))),data=d,weights=clone_weight,family=quasibinomial(),control=glm.control(maxit=100))
   e<-predict(f,type='response');ow<-d$clone_weight*ifelse(d$high==1,1-e,e)
   x<-model.matrix(f);x<-x[,colnames(x)!='(Intercept)',drop=FALSE]
   for(j in seq_len(ncol(x))){
     sdpooled<-sqrt((weighted_var(x[d$high==1,j],d$clone_weight[d$high==1])+weighted_var(x[d$high==0,j],d$clone_weight[d$high==0]))/2)
     if(!is.finite(sdpooled)||sdpooled<1e-12)next
     rawdiff<-weighted_average(x[d$high==1,j],d$clone_weight[d$high==1])-weighted_average(x[d$high==0,j],d$clone_weight[d$high==0])
     wdiff<-weighted_average(x[d$high==1,j],ow[d$high==1])-weighted_average(x[d$high==0,j],ow[d$high==0])
     balances[[length(balances)+1]]<-data.table(dataset=nm,specification=spec,covariate=colnames(x)[j],raw_smd=rawdiff/sdpooled,overlap_smd=wdiff/sdpooled)
   }
   cw<-data.table(clone=d$clone_id,high=d$high,w=ow)[,.(w=sum(w)),by=.(clone,high)]
   ess<-cw[,.(ess=sum(w)^2/sum(w^2),weight=sum(w)),by=high]
   ci<-cluster_ci(y,d$high,ow,d$clone_id)
   overlaps[[length(overlaps)+1]]<-data.table(dataset=nm,specification=spec,units=nrow(d),clones=uniqueN(d$clone_id),stress_cutoff=threshold,high_clones=uniqueN(d[high==1]$clone_id),low_clones=uniqueN(d[high==0]$clone_id),ess_high=ess[high==1]$ess,ess_low=ess[high==0]$ess,extreme_propensity_units=sum(e<.05|e>.95),raw_high_minus_low=weighted_average(y[d$high==1],d$clone_weight[d$high==1])-weighted_average(y[d$high==0],d$clone_weight[d$high==0]),overlap_high_minus_low=ci['estimate'],conditional_lower=ci['lower'],conditional_upper=ci['upper'],converged=f$converged,rank_design=f$rank,design_columns=ncol(model.matrix(f)))
   ps_rows[[length(ps_rows)+1]]<-d[,.(dataset=nm,specification=spec,clone_id,future_condition,high,propensity=e,overlap_weight=ow)]
 }
 feat<-setdiff(intersect(desired_features,names(d)),stress_axes)
 branch<-grep('^p_branch_',names(d),value=TRUE)
 full_nn<-neighbours(d,feat)
 # Record how selecting deeper clones changes condition and library composition.
 for(m in c(1,2,5,10)){
  ds<-d[descendant_cells>=m]
  sel<-ds[,.(units=.N,clones=uniqueN(clone_id),mean_stress=mean(stress_state),mean_log_detected_genes=mean(log_detected_genes),mean_mitochondrial_fraction=mean(mitochondrial_fraction)),by=.(future_condition,dominant_library)];sel[,`:=`(dataset=nm,minimum_descendants=m)];selection[[length(selection)+1]]<-sel
  # Estimand A: same original geometry, only the evaluated clones are selected.
  ss<-fit_one(ds,ds$local_effective_rank_stress_state_loo,TRUE)
  ss[,`:=`(dataset=nm,minimum_descendants=m,geometry='full_population_geometry_selected_queries')];subset_slopes[[length(subset_slopes)+1]]<-ss
  # Estimand B: selected clones define both queries and re-fitted geometry.
  if(m>1){nn<-neighbours(ds,feat);g<-geometry(as.matrix(ds[,..branch]),nn);ss<-fit_one(ds,g$rank,TRUE);ss[,`:=`(dataset=nm,minimum_descendants=m,geometry='selected_population_geometry_selected_queries')];subset_slopes[[length(subset_slopes)+1]]<-ss}
  if(m>=2){
   # Estimand C: preserve the original PCA and evaluate the original queries,
   # changing only the donor pool to deep clones; retain original k when possible.
   allrank<-numeric(nrow(d));allk<-integer(nrow(d))
   for(co in unique(d$future_condition)){
     rr<-which(d$future_condition==co); eligible<-which(d$future_condition==co & d$descendant_cells>=m)
     x<-scale(as.matrix(d[rr,..feat]));x[!is.finite(x)]<-0
     zz<-prcomp(x,rank.=min(15,ncol(x),nrow(x)-1),center=FALSE,scale.=FALSE)$x
     local_pool<-match(eligible,rr);y_pool<-as.matrix(d[eligible,..branch]);k0<-full_nn$k[rr[1]]
     for(i in seq_along(rr)){
       distance<-rowSums((sweep(zz[local_pool,,drop=FALSE],2,zz[i,],'-'))^2)
       distance[eligible==rr[i]]<-Inf
       k<-min(k0,sum(is.finite(distance)));idx<-order(distance)[seq_len(k)]
       cv<-cov(y_pool[idx,,drop=FALSE]);ev<-pmax(eigen(cv,symmetric=TRUE,only.values=TRUE)$values,0);pp<-ev/max(sum(ev),1e-12)
       allrank[rr[i]]<-if(sum(ev)>0)exp(-sum(pp[pp>0]*log(pp[pp>0])))else 0;allk[rr[i]]<-k
     }
   }
   ss<-fit_one(d,allrank,TRUE);ss[,`:=`(dataset=nm,minimum_descendants=m,geometry='fixed_full_PCA_original_queries_deep_donors',k_min=min(allk),k_max=max(allk))];donor_results[[length(donor_results)+1]]<-ss
  }
 }
 # Original full-geometry within-condition coefficients; no between-condition contrast.
 for(co in unique(d$future_condition)){
  ds<-d[future_condition==co];ss<-fit_one(ds,ds$local_effective_rank_stress_state_loo,TRUE);ss[,`:=`(dataset=nm,future_condition=co)];condition_slopes[[length(condition_slopes)+1]]<-ss
 }
}
bal<-rbindlist(balances);fwrite(bal,file.path(out,'quality_balance_diagnostics.tsv'),sep='\t')
ov<-rbindlist(overlaps);ov<-merge(ov,bal[,.(max_abs_raw_smd=max(abs(raw_smd)),max_abs_overlap_smd=max(abs(overlap_smd))),by=.(dataset,specification)]);fwrite(ov,file.path(out,'quality_overlap_summary.tsv'),sep='\t')
fwrite(rbindlist(ps_rows),file.path(out,'quality_overlap_weights.tsv'),sep='\t')
fwrite(rbindlist(selection),file.path(out,'descendant_selection_composition.tsv'),sep='\t')
fwrite(rbindlist(subset_slopes),file.path(out,'fixed_geometry_subset_slopes.tsv'),sep='\t')
fwrite(rbindlist(donor_results),file.path(out,'fixed_PCA_deep_donor_slopes.tsv'),sep='\t')
fwrite(rbindlist(condition_slopes),file.path(out,'original_geometry_condition_slopes.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(ov);print(rbindlist(subset_slopes));print(rbindlist(donor_results));print(rbindlist(condition_slopes))

# Refit overlap weights within each clone-cluster bootstrap; geometry stays fixed.
cache<-readRDS(file.path(pkg,'inputs','lineage_inputs.rds'));a<-list();warn<-list()
for(nm in names(cache$data)){
 d<-copy(cache$data[[nm]]);d[,clone_weight:=1/.N,by=clone_id];cut<-median(unique(d,by='clone_id')$stress_state);d[,high:=as.integer(stress_state>cut)]
 for(spec in c('quality_only','quality_and_other_states')){
  rhs<-c('factor(baseline_n)','splines::ns(log_descendant_n,df=2)','splines::ns(log_detected_genes,df=3)','splines::ns(mitochondrial_fraction,df=3)')
  if(spec=='quality_and_other_states')rhs<-c(rhs,'cell_cycle_state','chromatin_regulator_proxy','metabolic_state')
  if(uniqueN(d$future_condition)>1)rhs<-c(rhs,'factor(future_condition)')
  if(uniqueN(d$dominant_library)>1)rhs<-c(rhs,'factor(dominant_library)')
  x<-model.matrix(as.formula(paste('~',paste(rhs,collapse='+'))),d)
  ids<-unique(d$clone_id);indices<-split(seq_len(nrow(d)),d$clone_id)
  set.seed(202611033+100*match(nm,names(cache$data))+match(spec,c('quality_only','quality_and_other_states')))
  for(b in seq_len(499)){
    z<-unlist(indices[as.character(sample(ids,length(ids),replace=TRUE))],use.names=FALSE)
    fit<-suppressWarnings(glm.fit(x[z,,drop=FALSE],d$high[z],weights=d$clone_weight[z],family=quasibinomial(),control=glm.control(maxit=100)))
    e<-fit$fitted.values;high<-d$high[z];w<-d$clone_weight[z]*ifelse(high==1,1-e,e);y<-d$local_effective_rank_stress_state_loo[z]
    diff<-weighted.mean(y[high==1],w[high==1])-weighted.mean(y[high==0],w[high==0])
    a[[length(a)+1]]<-data.table(dataset=nm,specification=spec,replicate=b,rank_difference=diff,converged=fit$converged)
  }
 }
}
a<-rbindlist(a);fwrite(a,file.path(out,'quality_overlap_bootstrap.tsv'),sep='\t')
s<-a[,.(replicates=.N,valid=sum(is.finite(rank_difference)&converged),bootstrap_lower=quantile(rank_difference,.025),bootstrap_upper=quantile(rank_difference,.975),bootstrap_median=median(rank_difference)),by=.(dataset,specification)]
fwrite(s,file.path(out,'quality_overlap_bootstrap_summary.tsv'),sep='\t');print(s)

# Paired-clone context comparisons in the same source experiment.
d<-readRDS(file.path(pkg,'inputs','lineage_inputs.rds'))$data$larry_cytokine
branches<-grep('^p_branch_',names(d),value=TRUE)
results<-list();records<-list();branchmeans<-list()
target<-c('M-csf'='p_branch_monocyte','G-csf'='p_branch_neutrophil','Epo'='p_branch_erythroid')
jsd<-function(p,q){m<-(p+q)/2;p0<-p>0;q0<-q>0;.5*(sum(p[p0]*log(p[p0]/m[p0]))+sum(q[q0]*log(q[q0]/m[q0])))}
for(co in names(target)){
 a<-d[future_condition=='full cocktail'];b<-d[future_condition==co]
 shared<-intersect(a$clone_id,b$clone_id);a<-a[match(shared,clone_id)];b<-b[match(shared,clone_id)]
 stopifnot(identical(a$clone_id,b$clone_id),max(abs(a$stress_state-b$stress_state))<1e-12,max(abs(a$mitochondrial_fraction-b$mitochondrial_fraction))<1e-12,max(abs(a$detected_genes-b$detected_genes))<1e-12)
 for(depth in c(1,2)){
  keep<-a$descendant_cells>=depth & b$descendant_cells>=depth;aa<-a[keep];bb<-b[keep]
  ya<-as.matrix(aa[,..branches]);yb<-as.matrix(bb[,..branches]);delta<-yb-ya;n<-nrow(delta)
  j<-match(target[[co]],branches);v<-delta[,j];obs<-sum(colMeans(delta)^2)
  set.seed(202611034+100*match(co,names(target))+depth)
  signs<-matrix(sample(c(-1,1),9999*n,replace=TRUE),nrow=9999,ncol=n)
  null<-rowSums((signs%*%delta/n)^2);p<-(1+sum(null>=obs-1e-12))/10000
  ci<-quantile(replicate(1999,mean(sample(v,length(v),replace=TRUE))),c(.025,.975))
  jj<-vapply(seq_len(n),function(i)jsd(ya[i,],yb[i,]),numeric(1))
  results[[length(results)+1]]<-data.table(condition=co,reference='full cocktail',minimum_descendants_each=depth,paired_clones=n,mean_jsd_nats=mean(jj),median_jsd_nats=median(jj),l2_mean_probability_shift=sqrt(obs),target_branch=sub('p_branch_','',target[[co]]),mean_target_probability_reference=mean(ya[,j]),mean_target_probability_condition=mean(yb[,j]),paired_target_shift=mean(v),conditional_lower=ci[1],conditional_upper=ci[2],paired_sign_swap_diagnostic_p=p,baseline_stress_max_difference=max(abs(aa$stress_state-bb$stress_state)))
  bm<-data.table(branch=sub('p_branch_','',branches),mean_reference=colMeans(ya),mean_condition=colMeans(yb),mean_delta=colMeans(delta));bm[,`:=`(condition=co,minimum_descendants_each=depth,paired_clones=n)];branchmeans[[length(branchmeans)+1]]<-bm
  rec<-data.table(clone_id=aa$clone_id,reference_descendants=aa$descendant_cells,condition_descendants=bb$descendant_cells,target_delta=v,jsd_nats=jj);rec[,`:=`(condition=co,minimum_descendants_each=depth)];records[[length(records)+1]]<-rec
 }
}
res<-rbindlist(results);res[,diagnostic_bh:=p.adjust(paired_sign_swap_diagnostic_p,'BH'),by=minimum_descendants_each]
fwrite(res,file.path(out,'paired_clone_context_summary.tsv'),sep='\t')
fwrite(rbindlist(records),file.path(out,'paired_clone_context_records.tsv'),sep='\t')
fwrite(rbindlist(branchmeans),file.path(out,'paired_clone_branch_changes.tsv'),sep='\t');print(res)

writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Additional lineage/quality/context checks completed.\n')

# One explicit primary table keeps bootstrap uncertainty separate from fixed-weight diagnostics.
pt<-fread(file.path(out,'quality_overlap_summary.tsv'));setnames(pt,c('conditional_lower','conditional_upper'),c('fixed_weight_lower','fixed_weight_upper'))
bs<-fread(file.path(out,'quality_overlap_bootstrap_summary.tsv'))
fwrite(merge(pt,bs,by=c('dataset','specification')),file.path(out,'quality_overlap_results.tsv'),sep='\t')


#!/usr/bin/env Rscript
# Reconstruct original feature-specific geometry, permutation inference and transfer.
# Uses frozen clone-level source summaries; does not rerun raw-cell scoring.
suppressPackageStartupMessages({library(data.table);library(BiocNeighbors)})
args<-commandArgs(trailingOnly=TRUE)
pkg<-if(length(args)&&dir.exists(file.path(args[1],"inputs")))args[1]else"."
out<-if(length(args)>1&&!grepl("^--",args[2]))args[2]else"results/lineage_primary"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache<-readRDS(file.path(pkg,"inputs","lineage_inputs.rds"))
cluster_robust_coeftable <- function(fit, cluster) {
    beta <- coef(fit)
    estimable <- names(beta)[is.finite(beta)]
    x <- model.matrix(fit)[, estimable, drop = FALSE]
    e <- residuals(fit)
    bread <- solve(crossprod(x))
    score <- rowsum(x * e, group = cluster, reorder = FALSE)
    g <- nrow(score)
    n <- nrow(x)
    p <- ncol(x)
    correction <- if (g > 1L && n > p) {
        (g/(g - 1)) * ((n - 1)/(n - p))
    }
    else {
        1
    }
    vcov_cluster <- correction * bread %*% crossprod(score) %*% bread
    se <- sqrt(pmax(diag(vcov_cluster), 0))
    statistic <- beta[estimable]/se
    data.table(term = estimable, estimate = unname(beta[estimable]), std_error = se, statistic, p_value = 2 * 
        pt(-abs(statistic), df = pmax(g - 1L, 1L)), cluster_count = g)
}
calculate_local_geometry <- function(rows, clone, feature_set, branch_cols) {
    x <- as.matrix(clone[rows, ..feature_set])
    x <- scale(x)
    x[!is.finite(x)] <- 0
    pc <- prcomp(x, rank. = min(15L, ncol(x), nrow(x) - 1L), center = FALSE, scale. = FALSE)
    z <- pc$x
    y <- as.matrix(clone[rows, ..branch_cols])
    y[!is.finite(y)] <- 0
    y <- y/pmax(rowSums(y), 1)
    k <- min(60L, max(15L, floor(sqrt(nrow(z)))), nrow(z) - 1L)
    if (k < 2L) 
        return(NULL)
    nn <- BiocNeighbors::findKNN(z, k = k, BNPARAM = KmknnParam())
    rbindlist(lapply(seq_len(nrow(z)), function(i) {
        idx <- nn$index[i, ]
        yy <- y[idx, , drop = FALSE]
        cov_y <- stats::cov(yy)
        eig <- pmax(eigen(cov_y, symmetric = TRUE, only.values = TRUE)$values, 0)
        p <- eig/pmax(sum(eig), 1e-12)
        effective_rank <- if (sum(eig) > 0) {
            exp(-sum(p[p > 0] * log(p[p > 0])))
        }
        else {
            0
        }
        local_prob <- colMeans(yy)
        target_energy <- 1/(diag(cov_y) + 1e-04)
        data.table(row_id = rows[i], neighbourhood_size = k, local_effective_rank = effective_rank, dominant_direction_fraction = max(p, 
            na.rm = TRUE), unreachable_fates = sum(local_prob < 0.05), local_transition_entropy = -sum(local_prob[local_prob > 
            0] * log(local_prob[local_prob > 0])), minimum_local_control_energy = min(target_energy, 
            na.rm = TRUE), maximum_local_control_energy = max(target_energy, na.rm = TRUE), local_reachability = max(local_prob, 
            na.rm = TRUE))
    }))
}
partial_clone_permutation_p <- function(model_dt, outcome, focal, state_features, condition_adjustment, 
    permutations = 999L) {
    other <- setdiff(state_features, focal)
    clone_cov <- unique(model_dt[, c("clone_id", state_features), with = FALSE], by = "clone_id")
    reduced <- lm(as.formula(paste(focal, "~", paste(other, collapse = " + "))), data = clone_cov)
    focal_fitted <- fitted(reduced)
    focal_residual <- residuals(reduced)
    full_formula <- as.formula(paste(outcome, "~", paste(state_features, collapse = " + "), condition_adjustment))
    observed <- coef(lm(full_formula, data = model_dt))[[focal]]
    null <- numeric(permutations)
    for (b in seq_len(permutations)) {
        permuted <- focal_fitted + sample(focal_residual, replace = FALSE)
        perm_map <- setNames(permuted, clone_cov$clone_id)
        perm_dt <- copy(model_dt)
        perm_dt[, `:=`((focal), unname(perm_map[as.character(clone_id)]))]
        null[b] <- coef(lm(full_formula, data = perm_dt))[[focal]]
    }
    (1 + sum(abs(null) >= abs(observed), na.rm = TRUE))/(permutations + 1)
}
ridge_transfer <- function(train, test, outcome, lambda = 10) {
    xtr <- as.matrix(train[, ..common_features])
    xte <- as.matrix(test[, ..common_features])
    mu <- colMeans(xtr)
    sig <- apply(xtr, 2, sd)
    sig[!is.finite(sig) | sig < 1e-08] <- 1
    xtr <- sweep(sweep(xtr, 2, mu, "-"), 2, sig, "/")
    xte <- sweep(sweep(xte, 2, mu, "-"), 2, sig, "/")
    xtr[!is.finite(xtr)] <- 0
    xte[!is.finite(xte)] <- 0
    y <- train[[outcome]]
    beta <- solve(crossprod(cbind(1, xtr)) + diag(c(0, rep(lambda, ncol(xtr)))), crossprod(cbind(1, xtr), 
        y))
    as.numeric(cbind(1, xte) %*% beta)
}

state_features<-c('cell_cycle_state','stress_state','chromatin_regulator_proxy','metabolic_state','mitochondrial_fraction')
excluded_axis_features<-list(
 cell_cycle_state=c('hallmark__E2F_TARGETS','hallmark__G2M_CHECKPOINT'),
 stress_state=c('hallmark__HYPOXIA','hallmark__REACTIVE_OXYGEN_SPECIES_PATHWAY','hallmark__UNFOLDED_PROTEIN_RESPONSE'),
 chromatin_regulator_proxy=c('collectri__Ezh2','collectri__Dnmt1','collectri__Dnmt3a','collectri__Hdac1','collectri__Kdm5a'),
 metabolic_state='hallmark__OXIDATIVE_PHOSPHORYLATION',mitochondrial_fraction=character())
objects<-lapply(cache$data,function(d){
 d<-copy(d)
 if(!'transition_entropy_primary'%in%names(d))stop('Frozen input requires the source-precision transition_entropy_primary column for exact entropy-rank transfer reproduction.')
 # The legacy entropy column came through a 15-digit display TSV and collapsed
 # a few entropy ties. This column is frozen directly from the source fate table.
 d[,transition_entropy:=transition_entropy_primary]
 d
})
common_features<-Reduce(intersect,lapply(objects,function(d)intersect(cache$features,names(d))))
associations<-list();geometries<-list();agreement<-list()
set.seed(1905)
for(nm in names(objects)){
 clone<-objects[[nm]];setorder(clone,clone_id,future_condition);clone[,row_id:=.I]
 features<-intersect(cache$features,names(clone));branch_cols<-grep('^p_branch_',names(clone),value=TRUE)
 condition_rows<-split(clone$row_id,clone$future_condition)
 # Reconstruct the original full geometry before its feature-specific variants.
 full<-rbindlist(lapply(condition_rows,calculate_local_geometry,clone=clone,feature_set=features,branch_cols=branch_cols))
 full<-full[match(clone$row_id,row_id)]
 agreement[[length(agreement)+1]]<-data.table(dataset_id=nm,quantity='full_effective_rank',maximum_absolute_difference=max(abs(full$local_effective_rank-clone$local_effective_rank)))
 for(focal in state_features){
   local<-rbindlist(lapply(condition_rows,calculate_local_geometry,clone=clone,feature_set=setdiff(features,excluded_axis_features[[focal]]),branch_cols=branch_cols))
   local<-local[match(clone$row_id,row_id)]
   if(focal=='stress_state'){
    agreement[[length(agreement)+1]]<-data.table(dataset_id=nm,quantity='stress_excluded_effective_rank',maximum_absolute_difference=max(abs(local$local_effective_rank-clone$local_effective_rank_stress_state_loo)))
    agreement[[length(agreement)+1]]<-data.table(dataset_id=nm,quantity='stress_excluded_inverse_variance',maximum_absolute_difference=max(abs(local$minimum_local_control_energy-clone$minimum_local_control_energy_stress_state_loo)))
   }
   gg<-local[,.(row_id,local_effective_rank,minimum_local_control_energy)];gg[,`:=`(dataset_id=nm,clone_id=clone$clone_id,future_condition=clone$future_condition,excluded_axis=focal)];geometries[[length(geometries)+1]]<-gg
   for(outcome_base in c('local_effective_rank','minimum_local_control_energy'))clone[,(paste0(outcome_base,'_',focal,'_loo')):=local[[outcome_base]]]
 }
 for(outcome_base in c('local_effective_rank','minimum_local_control_energy')){
  for(focal in state_features){
   outcome<-paste0(outcome_base,'_',focal,'_loo')
   model_dt<-clone[complete.cases(clone[,c(outcome,state_features),with=FALSE])]
   condition_adjustment<-if(uniqueN(model_dt$future_condition)>1L)' + factor(future_condition)'else''
   form<-as.formula(paste(outcome,'~',paste(state_features,collapse=' + '),condition_adjustment))
   fit<-lm(form,data=model_dt);co<-cluster_robust_coeftable(fit,model_dt$clone_id);focal_row<-co[term==focal]
   permutation_seed<-1905+match(nm,names(objects))*100+match(focal,state_features)*10+match(outcome_base,c('local_effective_rank','minimum_local_control_energy'))
   set.seed(permutation_seed)
   permutation_p<-partial_clone_permutation_p(model_dt,outcome,focal,state_features,condition_adjustment,permutations=999L)
   focal_row[,`:=`(dataset_id=nm,outcome=outcome_base,neighborhood_excluded_axis=focal,asymptotic_p_value=p_value,p_value=permutation_p,inference_method='partial_clone_permutation_with_axis_excluded_neighborhood',permutations=999L,n=nrow(model_dt),permutation_seed=permutation_seed)]
   associations[[length(associations)+1]]<-focal_row
   cat(nm,outcome_base,focal,'P=',permutation_p,'\n')
  }
 }
 objects[[nm]]<-clone
}
a<-rbindlist(associations);a[,p_adj_bh:=p.adjust(p_value,'BH'),by=.(dataset_id,outcome)]
fwrite(a,file.path(out,'state_associations.tsv'),sep='\t')
fwrite(rbindlist(geometries),file.path(out,'primary_feature_geometry.tsv'),sep='\t')
ag<-rbindlist(agreement);fwrite(ag,file.path(out,'primary_geometry_reconstruction.tsv'),sep='\t');stopifnot(all(ag$maximum_absolute_difference<1e-9))

transfer<-list()
for(train_nm in names(objects)){
 test_nm<-setdiff(names(objects),train_nm);train<-objects[[train_nm]];test<-objects[[test_nm]]
 for(outcome in c('p_differentiation','transition_entropy','p_rare_fate')){
  pred<-ridge_transfer(train,test,outcome);obs<-test[[outcome]]
  transfer[[length(transfer)+1]]<-data.table(trained_on=train_nm,held_out=test_nm,outcome=outcome,n_test=length(obs),spearman=cor(obs,pred,method='spearman',use='complete.obs'),mse=mean((obs-pred)^2,na.rm=TRUE))
 }
}
tr<-rbindlist(transfer);fwrite(tr,file.path(out,'leave_dataset_transfer.tsv'),sep='\t')
aux<-readRDS(file.path(pkg,'inputs','auxiliary.rds'))
compare<-merge(tr,aux$larry_transfer,by=c('trained_on','held_out','outcome'))
transfer_error<-data.table(quantity=c('transfer_Spearman','transfer_MSE','transfer_test_n'),maximum_absolute_difference=c(max(abs(compare$spearman.x-compare$spearman.y)),max(abs(compare$mse.x-compare$mse.y)),max(abs(compare$n_test.x-compare$n_test.y))))
fwrite(transfer_error,file.path(out,'transfer_reconstruction.tsv'),sep='\t');stopifnot(nrow(compare)==6L,all(transfer_error$maximum_absolute_difference<1e-9))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(a);print(tr);cat('Original five-state/two-outcome permutation inference and transfer reconstructed.\n')
