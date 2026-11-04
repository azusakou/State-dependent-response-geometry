#!/usr/bin/env Rscript
# Derive compact count inputs directly from the deposited MultiPerturb object.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SingleCellExperiment)})
derive_primary <- function(sce) {
 safe_cosine<-function(x,y,eps=1e-10){den<-sqrt(sum(x*x))*sqrt(sum(y*y));if(!is.finite(den)||den<eps)return(NA_real_);sum(x*y)/den}
guide_id <- as.character(colData(sce)$guide_id)
guide_gene <- as.character(colData(sce)$guide_gene)
context <- sub("_.*$", "", guide_id)
context <- fifelse(context == "ATRT", "ATRT_human", "MNT_mouse")
guide_levels <- unique(guide_id)
guide_code <- match(guide_id, guide_levels)
n_guide <- length(guide_levels)
design <- sparseMatrix(
  i = seq_along(guide_code), j = guide_code, x = 1,
  dims = c(length(guide_code), n_guide)
)

guide_meta <- unique(data.table(
  guide_id = guide_id,
  guide_gene = guide_gene,
  context = context
), by = "guide_id")
guide_meta <- guide_meta[match(guide_levels, guide_id)]
guide_meta[, n_cells := tabulate(guide_code, nbins = n_guide)]
guide_meta[, is_control := grepl("^non-targeting", guide_gene)]

rna <- assay(sce, "counts")
rna_cell_library <- Matrix::colSums(rna)
rna_sum <- rna %*% design
rna_library <- as.numeric(crossprod(design, rna_cell_library))
rna_norm <- as.matrix(rna_sum)
rna_norm <- log1p(sweep(rna_norm, 2, pmax(rna_library, 1), "/") * 1e6)
rna_gene <- as.character(rowData(sce)$gene_symbol)
rownames(rna_norm) <- make.unique(ifelse(
  is.na(rna_gene) | !nzchar(rna_gene),
  as.character(rowData(sce)$ensembl_gene),
  rna_gene
))

atac_sce <- altExp(sce, "ATAC")
atac <- assay(atac_sce, "counts")
atac_cell_library <- Matrix::colSums(atac)
peak_gene <- as.character(rowData(atac_sce)$gene_symbol)
valid_peak <- !is.na(peak_gene) & nzchar(peak_gene) & peak_gene != "NA"
peak_gene_levels <- unique(peak_gene[valid_peak])
peak_map <- sparseMatrix(
  i = match(peak_gene[valid_peak], peak_gene_levels),
  j = which(valid_peak), x = 1,
  dims = c(length(peak_gene_levels), nrow(atac))
)
atac_gene_sum <- peak_map %*% (atac %*% design)
atac_library <- as.numeric(crossprod(design, atac_cell_library))
atac_norm <- as.matrix(atac_gene_sum)
atac_norm <- log1p(sweep(atac_norm, 2, pmax(atac_library, 1), "/") * 1e6)
rownames(atac_norm) <- peak_gene_levels

response_against_context_control <- function(mat) {
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
                dimnames = dimnames(mat))
  for (ctx in unique(guide_meta$context)) {
    cols <- which(guide_meta$context == ctx)
    ctrl <- cols[guide_meta$is_control[cols]]
    ctrl_mean <- rowSums(mat[, ctrl, drop = FALSE] *
      rep(guide_meta$n_cells[ctrl], each = nrow(mat))) /
      sum(guide_meta$n_cells[ctrl])
    out[, cols] <- sweep(mat[, cols, drop = FALSE], 1, ctrl_mean, "-")
  }
  out
}

rna_delta <- response_against_context_control(rna_norm)
atac_delta <- response_against_context_control(atac_norm)

rna_var <- apply(rna_delta[, !guide_meta$is_control, drop = FALSE], 1, var,
                 na.rm = TRUE)
atac_var <- apply(atac_delta[, !guide_meta$is_control, drop = FALSE], 1, var,
                  na.rm = TRUE)
rna_keep <- names(sort(rna_var, decreasing = TRUE))[seq_len(min(2000L, sum(
  is.finite(rna_var)
)))]
atac_keep <- names(sort(atac_var, decreasing = TRUE))[seq_len(min(2000L, sum(
  is.finite(atac_var)
)))]

rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))
rna_movement <- apply(rna_delta[rna_keep, , drop = FALSE], 2, rms)
atac_movement <- apply(atac_delta[atac_keep, , drop = FALSE], 2, rms)

common_gene <- intersect(rownames(rna_delta), rownames(atac_delta))
common_gene <- common_gene[
  common_gene %in% union(rna_keep, atac_keep)
]
cross_modality_cosine <- vapply(seq_len(n_guide), function(j) {
  safe_cosine(rna_delta[common_gene, j], atac_delta[common_gene, j])
}, numeric(1))

guide_summary <- copy(guide_meta)
guide_summary[, `:=`(
  rna_movement = rna_movement,
  atac_movement = atac_movement,
  cross_modality_cosine = cross_modality_cosine,
  modality_log2_ratio = log2((rna_movement + 1e-6) /
                               (atac_movement + 1e-6)),
  mean_rna_library = rna_library / pmax(n_cells, 1),
  mean_atac_library = atac_library / pmax(n_cells, 1)
)]
for (ctx in unique(guide_meta$context)) {
  cell_idx <- which(context == ctx)
  control_cells <- cell_idx[grepl("^non-targeting", guide_gene[cell_idx])]
  guide_idx <- which(guide_meta$context == ctx)
  rna_mu <- mean(rna_cell_library[control_cells], na.rm = TRUE)
  rna_sd <- sd(rna_cell_library[control_cells], na.rm = TRUE)
  atac_mu <- mean(atac_cell_library[control_cells], na.rm = TRUE)
  atac_sd <- sd(atac_cell_library[control_cells], na.rm = TRUE)
  if (!is.finite(rna_sd) || rna_sd < 1e-8) rna_sd <- 1
  if (!is.finite(atac_sd) || atac_sd < 1e-8) atac_sd <- 1
  guide_summary[guide_idx, `:=`(
    nuisance_rna_library_z = (mean_rna_library - rna_mu) / rna_sd,
    nuisance_atac_library_z = (mean_atac_library - atac_mu) / atac_sd
  )]
}
guide_summary[, nuisance_separation_check :=
  is.finite(nuisance_rna_library_z) &
  is.finite(nuisance_atac_library_z) &
  pmax(abs(nuisance_rna_library_z), abs(nuisance_atac_library_z)) < 2
]
guide_summary[, eligible_guide := !is_control & n_cells >= 30L]


 list(rna_response=rna_delta[rna_keep,,drop=FALSE],atac_response=atac_delta[atac_keep,,drop=FALSE],guide_summary=guide_summary)
}

a<-commandArgs(TRUE);if(length(a)<3)stop('Usage: Rscript extract_multiome_counts.R source.RDS.gz AUTO_or_frozen_multiome.rds output.rds')
con<-gzfile(a[1],'rb');sce<-readRDS(gzcon(con));old<-if(toupper(a[2])=="AUTO")derive_primary(sce)else readRDS(a[2]);g<-as.data.table(old$guide_summary)
rna<-assay(sce,'counts');atac_sce<-altExp(sce,'ATAC');atac<-assay(atac_sce,'counts')
gi<-match(as.character(colData(sce)$guide_id),g$guide_id);stopifnot(!anyNA(gi))
rl<-Matrix::colSums(rna);al<-Matrix::colSums(atac)
rgenes<-as.character(rowData(sce)$gene_symbol);rnames<-make.unique(ifelse(is.na(rgenes)|!nzchar(rgenes),as.character(rowData(sce)$ensembl_gene),rgenes))
rkeep<-match(rownames(old$rna_response),rnames);stopifnot(!anyNA(rkeep))
pgenes<-as.character(rowData(atac_sce)$gene_symbol);valid<-!is.na(pgenes)&nzchar(pgenes)&pgenes!='NA';ag<-rownames(old$atac_response)
map_gene<-function(keep){jj<-which(keep&pgenes%in%ag);sparseMatrix(i=match(pgenes[jj],ag),j=jj,x=1,dims=c(length(ag),nrow(atac)))}
map_all<-map_gene(valid);ac<-map_all%*%atac
# The primary count/guide/quality gate precedes the descriptive cosine filter.
targets<-g[is_control==FALSE,.(n_guides=.N,n_cells=sum(n_cells),minimum_guide_cells=min(n_cells),all_nuisance=all(nuisance_separation_check)),by=.(context,guide_gene)]
targets[,count_quality_gate:=n_guides>=2&n_cells>=100&minimum_guide_cells>=15&all_nuisance]
selected_guides<-which(g$context=='ATRT_human'&(g$is_control|g$guide_gene%in%targets[count_quality_gate==TRUE,guide_gene]))
cc<-which(gi%in%selected_guides);ids<-colnames(sce);if(is.null(ids))ids<-sprintf('cell_%06d',seq_len(ncol(sce)))
# All target-specific RNA counts are used for a transcript-depletion proxy.
tg<-sort(unique(g$guide_gene[!g$is_control]));rr<-which(rgenes%in%tg)
trmap<-sparseMatrix(i=match(rgenes[rr],tg),j=rr,x=1,dims=c(length(tg),nrow(rna)))
tc<-trmap%*%rna
on_count<-numeric(ncol(rna));target_idx<-match(g$guide_gene[gi],tg);ok<-which(!is.na(target_idx));on_count[ok]<-as.numeric(tc[cbind(target_idx[ok],ok)])
nt<-which(g$context[gi]=='ATRT_human'&g$is_control[gi]);high<-which(g$context[gi]=='ATRT_human'&!g$is_control[gi]&g$n_cells[gi]>=15)
ot_cells<-sort(unique(c(nt,high)))
# Alternative deposited-link restrictions retain only genes with available peaks.
dist_tss<-suppressWarnings(as.numeric(rowData(atac_sce)$distanceToTSS))
des<-sparseMatrix(i=seq_len(ncol(sce)),j=gi,x=1,dims=c(ncol(sce),nrow(g)))
context_delta<-function(counts){s<-as.matrix(counts%*%des);lib<-as.numeric(crossprod(des,al));v<-log1p(sweep(s,2,pmax(lib,1),'/')*1e6);out<-v
 for(ctx in unique(g$context)){jj<-which(g$context==ctx);nn<-jj[g$is_control[jj]];mu<-as.numeric(v[,nn,drop=FALSE]%*%g$n_cells[nn]/sum(g$n_cells[nn]));out[,jj]<-sweep(v[,jj,drop=FALSE],1,mu,'-')};out}
windows<-list();peak_coverage<-list()
for(w in c(2000,10000)){mp<-map_gene(valid&is.finite(dist_tss)&abs(dist_tss)<=w);pc<-Matrix::rowSums(mp);delta<-context_delta(mp%*%atac);ii<-which(pc>0);rownames(delta)<-ag;windows[[as.character(w)]]<-delta[ii,,drop=FALSE];peak_coverage[[as.character(w)]]<-data.table(gene=ag,n_peaks=pc,tss_window=w)}
x<-list(rna_counts=rna[rkeep,cc,drop=FALSE],atac_counts=ac[,cc,drop=FALSE],rna_library=rl[cc],atac_library=al[cc],cell_guide=gi[cc],cell_id=ids[cc],guide=g,targets=targets,rna_features=rownames(old$rna_response),atac_features=ag,atac_tss_responses=windows,atac_peak_coverage=rbindlist(peak_coverage),on_target_cell_metadata=data.table(cell_id=ids[ot_cells],guide_index=gi[ot_cells],rna_library=rl[ot_cells],atac_library=al[ot_cells],target_rna_count=on_count[ot_cells]),control_target_counts=tc[,nt,drop=FALSE],control_rna_library=rl[nt],target_genes=tg,target_gene_feature_count=tabulate(match(rgenes,tg),nbins=length(tg)),control_cell_ids=ids[nt],source='GSE277747/GSM8528725_MPS.RDS.gz',seed=20261104)
rownames(x$rna_counts)<-x$rna_features;rownames(x$atac_counts)<-x$atac_features
saveRDS(x,a[3],compress='xz');if(length(a)>=4)saveRDS(old,a[4],compress='xz');cat('Stored',length(cc),'count-profile cells and',length(ot_cells),'on-target proxy cells; bytes',file.info(a[3])$size,'\n');print(targets[count_quality_gate==TRUE]);print(g[seq_len(nrow(g))%in%selected_guides,.(guides=.N,cells=sum(n_cells),min_rna=min(mean_rna_library),min_atac=min(mean_atac_library))]);print(rbindlist(peak_coverage)[,.(covered_genes=sum(n_peaks>0),peaks=sum(n_peaks)),by=tss_window])
