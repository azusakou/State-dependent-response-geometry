#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(TRUE)
input<-"inputs/multiome.rds"
out_dir<-"results/multiperturb_multiome"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(input)
rna_delta<-x$rna_response; atac_delta<-x$atac_response; guide_summary<-as.data.table(x$guide_summary)
rna_keep<-rownames(rna_delta); atac_keep<-rownames(atac_delta)
rms<-function(x) sqrt(mean(x^2,na.rm=TRUE))
safe_cosine<-function(x,y,eps=1e-10) {den<-sqrt(sum(x*x))*sqrt(sum(y*y));if(!is.finite(den)||den<eps) return(NA_real_);sum(x*y)/den}
target_repro <- guide_summary[is_control == FALSE, {
  cols <- .I
  rna_mean <- rowMeans(rna_delta[rna_keep, cols, drop = FALSE])
  atac_mean <- rowMeans(atac_delta[atac_keep, cols, drop = FALSE])
  rna_cos <- vapply(cols, function(j) {
    other <- setdiff(cols, j)
    if (!length(other)) return(NA_real_)
    reference <- rowMeans(
      rna_delta[rna_keep, other, drop = FALSE]
    )
    safe_cosine(rna_delta[rna_keep, j], reference)
  }, numeric(1))
  atac_cos <- vapply(cols, function(j) {
    other <- setdiff(cols, j)
    if (!length(other)) return(NA_real_)
    reference <- rowMeans(
      atac_delta[atac_keep, other, drop = FALSE]
    )
    safe_cosine(atac_delta[atac_keep, j], reference)
  }, numeric(1))
  list(
    n_guides = .N,
    n_cells = sum(n_cells),
    median_rna_guide_cosine = median(rna_cos, na.rm = TRUE),
    median_atac_guide_cosine = median(atac_cos, na.rm = TRUE),
    reproducibility_method = "leave_one_guide_out_cosine",
    rna_movement = rms(rna_mean),
    atac_movement = rms(atac_mean),
    source_control_cells = sum(
      guide_summary$n_cells[
        guide_summary$context == .BY$context & guide_summary$is_control
      ]
    ),
    nuisance_separation_check = all(nuisance_separation_check),
    passes_task_qc = sum(n_cells) >= 100L && .N >= 2L &&
      min(n_cells) >= 15L &&
      median(rna_cos, na.rm = TRUE) > 0.1 &&
      median(atac_cos, na.rm = TRUE) > 0.1 &&
      all(nuisance_separation_check)
  )
}, by = .(context, guide_gene)]

# Control roles are target-level labels and require independent-guide
# reproducibility in both modalities. They are not assigned to single guides.
target_roles <- copy(target_repro)
target_roles[, rna_high := rna_movement >= quantile(
  rna_movement[passes_task_qc], 0.75, na.rm = TRUE
), by = context]
target_roles[, atac_high := atac_movement >= quantile(
  atac_movement[passes_task_qc], 0.75, na.rm = TRUE
), by = context]
target_roles[, control_role := fifelse(
  !passes_task_qc, "insufficient_support",
  fifelse(rna_high & !atac_high, "state_mover",
    fifelse(!rna_high & atac_high, "barrier_breaker_candidate",
      fifelse(rna_high & atac_high, "coupled_mover", "weak_or_contextual")
    )
  )
)]


# Threshold sensitivity uses the same guide vectors and quality gates.
fwrite(guide_summary, file.path(out_dir, "multiome_guide_support.tsv"), sep="\t")
min_cells <- guide_summary[is_control == FALSE, .(min_guide_cells=min(n_cells)), by=.(context,guide_gene)]
audit <- merge(target_repro,min_cells,by=c("context","guide_gene"))
audit[, count_guide_quality_gate := n_cells>=100 & n_guides>=2 & min_guide_cells>=15 & nuisance_separation_check]
sens <- rbindlist(lapply(c(0,0.1,0.2,0.3,0.5), function(cut) {
 x <- audit[count_guide_quality_gate & median_rna_guide_cosine>cut & median_atac_guide_cosine>cut]
 data.table(cosine_threshold=cut, all_tasks=nrow(audit), count_quality_supported=sum(audit$count_guide_quality_gate), eligible_tasks=nrow(x), eligible_targets=uniqueN(x$guide_gene), ATRT_tasks=sum(x$context=="ATRT_human"), MNT_tasks=sum(x$context=="MNT_mouse"))
}))
fwrite(sens,file.path(out_dir,"multiome_threshold_sensitivity.tsv"),sep="\t")
# Movement cutoffs are descriptive quantiles recalculated among primary eligible tasks.
role_sens <- rbindlist(lapply(c(.5,.67,.75,.9),function(q) {
 x <- copy(audit[passes_task_qc == TRUE])
 x[, `:=`(rh=rna_movement>=quantile(rna_movement,q),ah=atac_movement>=quantile(atac_movement,q)),by=context]
 x[, role:=fifelse(rh&ah,"coupled_mover",fifelse(rh,"rna_dominant",fifelse(ah,"atac_dominant","weak_or_contextual")))]
 x[, .(quantile_cutoff=q,n_tasks=.N),by=role]
}))
fwrite(role_sens,file.path(out_dir,"multiome_role_sensitivity.tsv"),sep="\t")
# Gene-label randomization null uses distinct non-control guides within context,
# maintaining each tested group's guide count. It is not a negative-control null.
fast_logo <- function(gram,cols) {
 s <- gram[cols,cols,drop=FALSE]; n <- length(cols)
 if(n<2) return(NA_real_)
 num <- rowSums(s)-diag(s)
 den <- sqrt(pmax(diag(s),0)*pmax(sum(s)-2*rowSums(s)+diag(s),0))
 median(ifelse(den>1e-12,num/den,NA_real_),na.rm=TRUE)
}
G_r <- crossprod(rna_delta[rna_keep,,drop=FALSE]); G_a <- crossprod(atac_delta[atac_keep,,drop=FALSE])
set.seed(191103)
nulls <- rbindlist(lapply(seq_len(nrow(audit)),function(i) {
 d <- audit[i]; if(!d$count_guide_quality_gate) return(NULL)
 pool <- which(guide_summary$context==d$context & !guide_summary$is_control & guide_summary$n_cells>=15 & guide_summary$nuisance_separation_check)
 if(length(pool)<d$n_guides) return(NULL)
 sims <- replicate(999,{ cc<-sample(pool,d$n_guides,replace=FALSE); c(fast_logo(G_r,cc),fast_logo(G_a,cc))})
 data.table(context=d$context,guide_gene=d$guide_gene,n_guides=d$n_guides, label_null_joint_pass_0.1=mean(sims[1,]>.1&sims[2,]>.1), label_null_rna_95=quantile(sims[1,],.95), label_null_atac_95=quantile(sims[2,],.95),label_null_joint_p=(1+sum(pmin(sims[1,],sims[2,])>=min(d$median_rna_guide_cosine,d$median_atac_guide_cosine)))/1000)
}))
nulls[,label_null_joint_fdr:=p.adjust(label_null_joint_p,"BH")]
audit <- merge(audit,nulls,by=c("context","guide_gene","n_guides"),all.x=TRUE)
# Guide-resampling intervals describe construct variability only.
ints <- rbindlist(lapply(which(audit$passes_task_qc),function(i) {
 d <- audit[i]; cc <- which(guide_summary$context==d$context & guide_summary$guide_gene==d$guide_gene)
 sims <- replicate(999,{ c2<-sample(cc,length(cc),replace=TRUE); c(rms(rowMeans(rna_delta[rna_keep,c2,drop=FALSE])),rms(rowMeans(atac_delta[atac_keep,c2,drop=FALSE])))})
 data.table(context=d$context,guide_gene=d$guide_gene,rna_guide_boot_lo=quantile(sims[1,],.025),rna_guide_boot_hi=quantile(sims[1,],.975),atac_guide_boot_lo=quantile(sims[2,],.025),atac_guide_boot_hi=quantile(sims[2,],.975))
}))
audit <- merge(audit,ints,by=c("context","guide_gene"),all.x=TRUE)
fwrite(audit,file.path(out_dir,"multiome_target_audit.tsv"),sep="\t")
# Guide-level correlation is a descriptive depth diagnostic, not causal correction.
cors <- rbindlist(lapply(unique(guide_summary$context),function(ctx){
 x<-guide_summary[context==ctx & !is_control & n_cells>=15]
 if(nrow(x)<3) return(NULL)
 rbindlist(lapply(c("rna_movement","atac_movement"),function(y) rbindlist(lapply(c("n_cells","mean_rna_library","mean_atac_library"),function(z){data.table(context=ctx,response=y,covariate=z,n=nrow(x),spearman=cor(x[[y]],x[[z]],method="spearman",use="pairwise.complete.obs"))}))))
}))
fwrite(cors,file.path(out_dir,"multiome_depth_diagnostics.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"multiome_sessionInfo.txt"))
print(sens); print(role_sens); print(audit[passes_task_qc == TRUE,.(guide_gene,n_guides,n_cells,min_guide_cells,median_rna_guide_cosine,median_atac_guide_cosine,label_null_joint_p,label_null_joint_fdr)])

fwrite(target_repro,file.path(out_dir,"target_reproducibility.tsv"),sep="\t")
fwrite(target_roles,file.path(out_dir,"target_control_roles.tsv"),sep="\t")

fwrite(guide_summary,file.path(out_dir,"guide_summary.tsv"),sep="\t")
# Sampling-depth diagnostics within the primary eligible subset.
xg<-guide_summary[guide_gene%in%audit[passes_task_qc==TRUE,guide_gene]&is_control==FALSE]
xt<-audit[passes_task_qc==TRUE]
depth<-rbindlist(list(rbindlist(lapply(c('rna_movement','atac_movement'),function(r)rbindlist(lapply(c('n_cells','mean_rna_library','mean_atac_library'),function(z)data.table(scope='primary_eligible_guides',response=r,covariate=z,n=nrow(xg),spearman=cor(xg[[r]],xg[[z]],method='spearman')))))),rbindlist(lapply(c('rna_movement','atac_movement'),function(r)data.table(scope='primary_eligible_targets',response=r,covariate='n_cells',n=nrow(xt),spearman=cor(xt[[r]],xt$n_cells,method='spearman'))))))
fwrite(depth,file.path(out_dir,'multiome_eligible_depth_diagnostics.tsv'),sep='\t')
