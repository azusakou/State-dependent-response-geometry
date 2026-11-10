#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);project<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(project,'.Rlib'),.libPaths()));source(file.path(out,'phase2_helpers.R'))
p<-readRDS(file.path(out,'source_preparation.rds'));spec<-p$spec;meta<-p$source_cells;genes<-p$allgenes
coeff<-make_coefficients(p$universe,spec);x<-suppressWarnings(H5ADMatrix(p$source_path))
aggregation<-sparseMatrix(i=match(p$raw_gene,genes),j=seq_along(p$raw_gene),x=1,dims=c(length(genes),length(p$raw_gene)))
idx_phase2<-match(p$universe,genes);idx_original<-match(p$original_universe,genes);idx_upr<-match(p$original_upr,genes)
is_analysis<-meta$gene%in%p$analysis_genes|meta$is_control
units<-meta[is_analysis,.(n_cells=.N),by=.(gene,gem_group)];units[,unit_index:=.I]
unit_key<-paste(units$gene,units$gem_group,sep='|');unit_code<-match(paste(meta$gene,meta$gem_group,sep='|'),unit_key)
sum_score<-matrix(0,nrow(units),length(spec),dimnames=list(NULL,names(spec)))
controls<-meta[is_control==TRUE];control_map<-match(meta$cell_id,controls$cell_id)
cp<-matrix(NA_real_,nrow(controls),length(spec),dimnames=list(controls$cell_id,names(spec)));cn<-cp
ct<-matrix(NA_real_,nrow(controls),length(p$target_genes),dimnames=list(controls$cell_id,p$target_genes));raw_ct<-ct
lib_phase2<-numeric(nrow(meta));lib_original<-lib_phase2;assigned_raw<-rep(NA_real_,nrow(meta));original_upr<-lib_phase2;raw_upr_parts<-list()
for(st in seq.int(1,nrow(meta),by=1500L)){
 ii<-st:min(st+1499L,nrow(meta));b<-as.matrix(aggregation%*%as.matrix(x[,meta$source_row[ii],drop=FALSE]))
 l2<-colSums(b[idx_phase2,,drop=FALSE]);lo<-colSums(b[idx_original,,drop=FALSE]);stopifnot(all(l2>0),all(lo>0))
 lib_phase2[ii]<-l2;lib_original[ii]<-lo
 ti<-match(meta$gene[ii],genes);ok<-!is.na(ti);assigned_raw[ii[ok]]<-b[cbind(ti[ok],which(ok))]
 ur<-t(b[idx_upr,,drop=FALSE]);raw_upr_parts[[length(raw_upr_parts)+1L]]<-Matrix(ur,sparse=TRUE)
 original_upr[ii]<-rowMeans(log1p(ur*1e4/lo))
 use<-which(is_analysis[ii]);logb<-log1p(sweep(b[idx_phase2,use,drop=FALSE],2,l2[use],'/')*1e4)
 z<-score_log_counts(logb,p$universe,spec,coeff,meta$gene[ii[use]])
 codes<-unit_code[ii[use]];ss<-rowsum(z$scores,codes,reorder=FALSE);jj<-as.integer(rownames(ss));sum_score[jj,]<-sum_score[jj,]+ss
 ci<-which(meta$is_control[ii]);if(length(ci)){
  local<-match(ci,use);dest<-control_map[ii[ci]];cp[dest,]<-z$positive[local,,drop=FALSE];cn[dest,]<-z$negative[local,,drop=FALSE]
  tj<-match(p$target_genes,genes);found<-!is.na(tj)
  raw_ct[dest,found]<-t(b[tj[found],ci,drop=FALSE]);ct[dest,found]<-log1p(raw_ct[dest,found,drop=FALSE]*1e4/l2[ci])
 }
 cat('Source phase2/raw-cache scan',max(ii),'/',nrow(meta),'\n');flush.console()
}
stopifnot(all(is.finite(cp)),all(is.finite(cn)),all(is.finite(sum_score)))
raw_upr<-do.call(rbind,raw_upr_parts);colnames(raw_upr)<-p$original_upr;rownames(raw_upr)<-meta$cell_id
source_cache<-list(unit_metadata=units,unit_score_sums=sum_score,control_cells=controls,control_positive=cp,control_negative=cn,control_target_log=ct,analysis_genes=p$analysis_genes,spec=spec,universe=p$universe)
saveRDS(source_cache,file.path(out,'source_phase2_moment_cache.rds'),compress='gzip')
cells<-copy(meta);cells[,`:=`(context='K562',common_library=lib_original,phase2_library=lib_phase2,original_UPR=original_upr,target_raw_count=assigned_raw,target_measured_raw=gene%in%genes,target_in_original_universe=gene%in%p$original_universe,target_in_phase2_universe=gene%in%p$universe)]
target_info<-data.table(gene=p$target_genes,measured_raw=p$target_genes%in%genes,in_original_universe=p$target_genes%in%p$original_universe,in_phase2_universe=p$target_genes%in%p$universe)
raw_cache<-list(cells=cells,raw_upr=raw_upr,control_cells=controls,raw_control_targets=raw_ct,target_info=target_info,universe=p$original_universe,phase2_universe=p$universe,target_genes=p$target_genes,fixed_external_genes=p$fixed_external_genes)
saveRDS(raw_cache,file.path(out,'original91_source_raw_cache.rds'),compress='gzip')
e<-gene_effects_from_cache(source_cache,spec);stopifnot(all(is.finite(e$gene$deltaZ)))
fwrite(e$gene,file.path(out,'source_phase2_gene_effects.tsv.gz'),sep='\t');fwrite(e$gem,file.path(out,'source_phase2_GEM_effects.tsv.gz'),sep='\t')
fwrite(target_info,file.path(out,'raw_target_availability.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'source_scoring_session_info.txt'))
cat('Source phase2 effects and exact original91/7848 raw cache complete; RPE1 X unread.\n')
