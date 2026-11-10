#!/usr/bin/env Rscript
# Outcome-blind code preparation; execution requires the complete source lock.
args<-commandArgs(TRUE);stopifnot(length(args)>=3L)
project<-normalizePath(args[1]);out<-normalizePath(args[2]);manifest_file<-normalizePath(args[3])
.libPaths(c(file.path(project,'.Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(hdf5r);library(digest)})
locks<-fread(manifest_file);stopifnot(all(c('file','sha256')%in%names(locks)))
paths<-ifelse(grepl('^/',locks$file),locks$file,file.path(out,locks$file))
stopifnot(all(file.exists(paths)),all(vapply(paths,digest,character(1),file=TRUE,algo='sha256')==locks$sha256))
script_file<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))
stopifnot(any(normalizePath(paths)==script_file))
source(file.path(out,'phase2_helpers.R'))
prep<-readRDS(file.path(out,'source_preparation.rds'));spec<-readRDS(file.path(out,'phase2_axis_spec.rds'));universe<-prep$universe;upr<-spec$UPR$positive;analysis_genes<-prep$analysis_genes;targets<-prep$target_genes
stopifnot(length(universe)==6780L,length(upr)==90L,length(analysis_genes)==438L,length(spec)==105L)
coeff<-make_coefficients(universe,spec)
meta_all<-fread(file.path(out,'RPE1_author_roster_metadata_only.tsv.gz'))
meta_all<-meta_all[is_control==TRUE|gene%in%analysis_genes];setorder(meta_all,source_row)
stopifnot(!anyDuplicated(meta_all$source_row),all(is.finite(meta_all$UMI_count)&meta_all$UMI_count>0))
featuremap<-fread(file.path(project,'paper19opp-revision2/reconstruction/external_resources/rpe1_metadata_only_feature_mapping.tsv'))
rawgenes<-featuremap$gene;allgenes<-unique(rawgenes)
raw_to_common<-which(rawgenes%in%universe)
agg<-sparseMatrix(i=match(rawgenes[raw_to_common],universe),j=seq_along(raw_to_common),x=1,dims=c(length(universe),length(raw_to_common)))
stopifnot(all(universe%in%rawgenes),sum(duplicated(rawgenes))==1L)
ctrl_index<-which(meta_all$is_control);n<-nrow(meta_all);nc<-length(ctrl_index);np<-length(spec)
positive<-negative<-scores<-matrix(NA_real_,n,np,dimnames=list(NULL,names(spec)))
raw_upr<-matrix(NA_real_,n,length(upr),dimnames=list(NULL,upr))
control_targets<-matrix(NA_real_,nc,length(targets),dimnames=list(NULL,targets))
control_rowmap<-integer(n);control_rowmap[ctrl_index]<-seq_along(ctrl_index)
for(nm in c('common_library','common_detected_genes','matrix_library','assigned_target_count','assigned_target_logCP10k'))meta_all[,(nm):=NA_real_]
source_path<-'/Volumes/wang/paper6yhn/raw/replogle_rpe1/rpe1_raw_singlecell_01.h5ad'
f<-hdf5r::H5File$new(source_path,'r')
stopifnot(identical(as.integer(f[['X']]$dims),c(length(rawgenes),247914L)))
# All lock checks above run before the first access to any RPE1 X values.
writeLines(c(paste('UTC',format(Sys.time(),tz='UTC',usetz=TRUE)),paste('verified_manifest',manifest_file),paste('manifest_sha256',digest(manifest_file,file=TRUE,algo='sha256'))),file.path(out,'RPE1_expression_read_started.txt'))
started<-Sys.time();cat('RPE1 authorized count scan started',format(started,tz='UTC'),'UTC\n');flush.console()
for(st in seq.int(1L,247914L,by=5000L)){
 en<-min(st+4999L,247914L);jj<-which(meta_all$source_row>=st&meta_all$source_row<=en);if(!length(jj))next
 block<-f[['X']][,st:en];within<-meta_all$source_row[jj]-st+1L
 b<-as.matrix(agg%*%block[raw_to_common,within,drop=FALSE]);stopifnot(all(is.finite(b)),all(b>=0),max(abs(b-round(b)))<1e-6)
 lib<-colSums(b);valid<-which(is.finite(lib)&lib>0)
 set(meta_all,i=jj,j='common_library',value=lib);set(meta_all,i=jj,j='common_detected_genes',value=colSums(b>0));set(meta_all,i=jj,j='matrix_library',value=colSums(block[,within,drop=FALSE]))
 raw_upr[jj,]<-t(b[match(upr,universe),,drop=FALSE])
 genej<-meta_all$gene[jj];self<-match(genej,universe);targetraw<-rep(NA_real_,length(jj));ok<-which(!is.na(self));targetraw[ok]<-b[cbind(self[ok],ok)]
 outside<-which(is.na(self)&!meta_all$is_control[jj]);if(length(outside))for(g in unique(genej[outside])){kk<-outside[genej[outside]==g];ix<-which(rawgenes==g);if(length(ix))targetraw[kk]<-colSums(block[ix,within[kk],drop=FALSE])}
 set(meta_all,i=jj,j='assigned_target_count',value=targetraw)
 if(length(valid)){
  lg<-log1p(sweep(b[,valid,drop=FALSE],2,lib[valid],'/')*1e4)
  sc<-score_log_counts(lg,universe,spec,coeff,genej[valid]);scores[jj[valid],]<-sc$scores;positive[jj[valid],]<-sc$positive;negative[jj[valid],]<-sc$negative
  set(meta_all,i=jj[valid],j='assigned_target_logCP10k',value=log1p(1e4*targetraw[valid]/lib[valid]))
 }
 cc<-which(meta_all$is_control[jj]);if(length(cc)){
  cri<-control_rowmap[jj[cc]]
  for(k in seq_along(targets)){
   g<-targets[k];mi<-match(g,universe)
   if(!is.na(mi))control_targets[cri,k]<-b[mi,cc]
   else{ix<-which(rawgenes==g);if(length(ix))control_targets[cri,k]<-colSums(block[ix,within[cc],drop=FALSE])}
  }
 }
 if(en%%50000L==0L){cat('RPE1 raw cells through',en,'retained/scored',sum(is.finite(meta_all$common_library)),'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'seconds\n');flush.console()}
}
f$close_all();valid<-is.finite(meta_all$common_library)&meta_all$common_library>0
fwrite(meta_all[!valid],file.path(out,'RPE1_unscorable_common_library_cells.tsv'),sep='\t')
stopifnot(all(is.finite(scores[valid,,drop=FALSE])),all(is.finite(raw_upr[valid,,drop=FALSE])))
# Gene roster remains the source-locked panel, even if a cell cannot be scored.
meta<-meta_all[valid];scores<-scores[valid,,drop=FALSE];positive<-positive[valid,,drop=FALSE];negative<-negative[valid,,drop=FALSE];raw_upr<-raw_upr[valid,,drop=FALSE]
ci<-which(meta$is_control);oldctrlvalid<-valid[ctrl_index];control_targets<-control_targets[oldctrlvalid,,drop=FALSE]
control_target_log<-log1p(sweep(control_targets,1,meta$common_library[ci],'/')*1e4)
units<-meta[is_control==FALSE,.(n_cells=.N),by=.(gene,gem_group)];setorder(units,gene,gem_group)
key<-paste(meta$gene,meta$gem_group,sep='|');uindex<-match(key,paste(units$gene,units$gem_group,sep='|'));treated<-which(!meta$is_control)
unit_sums<-rowsum(scores[treated,,drop=FALSE],uindex[treated],reorder=TRUE);stopifnot(nrow(unit_sums)==nrow(units))
cache<-list(study='RPE1',source_path=source_path,source_bytes=file.info(source_path)$size,universe=universe,spec=spec,analysis_genes=analysis_genes,target_genes=targets,
 cell_metadata=meta,cell_scores=scores,cell_UPR_raw_counts=raw_upr,control_cells=meta[ci],control_positive=positive[ci,,drop=FALSE],control_negative=negative[ci,,drop=FALSE],control_target_counts=control_targets,control_target_log=control_target_log,unit_metadata=units,unit_score_sums=unit_sums)
saveRDS(cache,file.path(out,'RPE1_phase2_cache.rds'),compress='gzip')
fwrite(meta,file.path(out,'RPE1_phase2_cell_metadata.tsv.gz'),sep='\t')
small<-cbind(meta,as.data.table(scores[,c('UPR','core_apoptosis','core_proliferation','core_stress','E2F'),drop=FALSE]));fwrite(small,file.path(out,'RPE1_phase2_key_cell_scores.tsv.gz'),sep='\t')
fwrite(meta[,.(scored_cells=.N,GEM_labels=uniqueN(gem_group),mean_common_library=mean(common_library),mean_common_detected_genes=mean(common_detected_genes),mean_recorded_UMI=mean(UMI_count),mean_mitopercent=mean(mitopercent)),by=.(gene,is_control)],file.path(out,'RPE1_phase2_measurement_accounting.tsv'),sep='\t')
writeLines(c(paste('source_locked_genes',length(analysis_genes)),paste('metadata_selected_cells',nrow(meta_all)),paste('scorable_cells',nrow(meta)),paste('unscorable_cells',sum(!valid)),paste('scorable_control_cells',length(ci)),paste('raw_feature_rows',length(rawgenes)),paste('unique_raw_symbols',length(allgenes)),paste('elapsed_seconds',as.numeric(difftime(Sys.time(),started,units='secs'))),'Fixed-denominator primary scores include assigned-target numerator exclusion; controls remain base per-cell scores until target-specific correction in the shared helper. Raw UPR and target counts permit exact real-UPR L-target sensitivity. No RPE1 outcome-dependent ranking, matching or eligibility modification was performed.'),file.path(out,'RPE1_extraction_summary.txt'))
cat('RPE1 extraction complete:',nrow(meta),'cells;',length(analysis_genes),'locked genes; elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n')
