#!/usr/bin/env Rscript
# Separately reimplement the published additive-direction statistics from
# preserved processed profiles. Usage: Rscript verify_norman.R BUNDLE OUT_DIR
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- readRDS(file.path(root, 'inputs', 'norman_summary.rds'))
layout <- attr(x, 'gem_storage', exact=TRUE)
dg <- array(NA_real_, layout$dim)
for (part in layout$parts) {
  file <- file.path(root, 'inputs', part$file)
  stopifnot(unname(tools::md5sum(file)) == part$md5)
  dg[,,part$indices] <- readRDS(file)
}
groups <- split(which(x$is_single & !is.na(x$single_gene)),
                x$single_gene[x$is_single & !is.na(x$single_gene)])
norm <- function(z) sqrt(sum(z*z))
single <- lapply(groups, function(ii)
  as.vector(x$delta[,ii,drop=FALSE] %*% x$n_cells[ii]/sum(x$n_cells[ii])))
single_gem <- lapply(groups, function(ii)
  sapply(seq_len(dim(dg)[3]), function(g) {
    w <- x$n_condition_gem[ii,g]
    as.vector(matrix(dg[,ii,g], nrow=dim(dg)[1]) %*% w/sum(w))
  }))
result <- do.call(rbind, lapply(which(x$is_pair), function(j) {
  genes <- strsplit(x$conditions[j], '_', fixed=TRUE)[[1]]
  stopifnot(length(genes) == 2, all(genes %in% names(single)))
  additive <- single[[genes[1]]] + single[[genes[2]]]
  u <- additive/max(norm(additive), 1e-10)
  residual <- x$delta[,j] - additive
  v <- colSums((dg[,j,] - single_gem[[genes[1]]] - single_gem[[genes[2]]]) * u)
  cosine <- function(z) sum(x$delta[,j]*z)/(norm(x$delta[,j])*norm(z))
  redundant <- max(cosine(single[[genes[1]]]),cosine(single[[genes[2]]])) > .8 &&
    norm(x$delta[,j]) <= 1.1 * max(norm(single[[genes[1]]]),norm(single[[genes[2]]]))
  test <- t.test(v)
  data.frame(combination=x$conditions[j], projection=sum(residual*u),
             gem_mean=mean(v), p=test$p.value,
             ci_low=test$conf.int[1], ci_high=test$conf.int[2],
             redundant=redundant)
}))
result$q <- p.adjust(result$p, 'BH')
result$class <- ifelse(result$q<.05 & result$gem_mean>0,'synergistic',
  ifelse(result$q<.05 & result$gem_mean<0,'antagonistic',
  ifelse(result$redundant,'redundant','context-dependent')))
reference <- read.delim(gzfile(file.path(root, 'tables', 'interactions.tsv.gz')))
reference <- reference[match(result$combination, reference$combination),]
differences <- c(projection=max(abs(result$projection-reference$interaction_projection)),
                 gem_mean=max(abs(result$gem_mean-reference$replicate_projection_mean)),
                 p=max(abs(result$p-reference$replicate_projection_p)),
                 q=max(abs(result$q-reference$replicate_projection_fdr)),
                 ci_low=max(abs(result$ci_low-reference$ci_low)),
                 ci_high=max(abs(result$ci_high-reference$ci_high)))
stopifnot(nrow(result)==131, all(result$class==reference$interaction_class),
          max(differences)<1e-10)
write.table(result,file.path(out,'norman_recomputed_interactions.tsv'),sep='\t',
            quote=FALSE,row.names=FALSE)
checks <- data.frame(metric=names(differences),max_absolute_discrepancy=as.vector(differences))
write.table(checks,file.path(out,'norman_verification.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(table(result$class));print(checks)
