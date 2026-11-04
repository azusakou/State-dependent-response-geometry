#!/usr/bin/env Rscript
# Optional read-only raw-count extraction for the independent train-fold RNA baseline.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(HDF5Array);library(DelayedMatrixStats);library(hdf5r)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<3)stop('Usage: Rscript extract_baseline_profiles.R BUNDLE_ROOT QC_H5AD_DIR OUTPUT_DIR [RESULTS_ROOT]')
bundle<-normalizePath(args[1]);raw<-normalizePath(args[2]);out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
result_root<-if(length(args)>=4)args[4]else file.path(bundle,'results')
source(file.path(bundle,'upstream/h5ad_annotation.R'))
d<-as.data.table(readRDS(file.path(bundle,'inputs/mixseq_modeling.rds')));reg<-fread(file.path(result_root,'task_registry/state_transition_tasks.tsv'));d<-d[pair_id%in%reg[dataset_id=='mixseq_2020'&status=='eligible',task_id]]
key<-function(ex,time,line)paste(ex,ifelse(is.na(time),'NA',as.character(time)),line,sep='|');d[,baseline_context:=key(experiment,time_hours,cell_line)]
files<-file.path(raw,c('DMSO_6hr_expt1.qc.h5ad','DMSO_24hr_expt1.qc.h5ad','DMSO_6hr_expt3.qc.h5ad','DMSO_24hr_expt3.qc.h5ad','DMSO_expt10.qc.h5ad'))
feature_lists<-lapply(files,function(f){h<-hdf5r::H5File$new(f,'r');g<-as.character(read_h5ad_var(h,'gene'));h$close_all();g});genes<-Reduce(intersect,feature_lists);stopifnot(!anyDuplicated(genes))
profiles<-list();metadata<-list();membership<-list()
for(k in seq_along(files)){
 f<-files[k];h<-hdf5r::H5File$new(f,'r');o<-data.table(cell_id=read_h5ad_obs(h,'cell_id'),cell_line=read_h5ad_obs(h,'cell_line'),experiment=as.integer(read_h5ad_obs(h,'experiment')),time_hours=as.numeric(read_h5ad_obs(h,'time_hours')),qc=as.logical(read_h5ad_obs(h,'passes_qc_final')),perturbation=read_h5ad_obs(h,'perturbation'));h$close_all()
 o[,baseline_context:=key(experiment,time_hours,cell_line)];o[,column_index:=.I];o<-o[qc==TRUE&tolower(perturbation)=='dmso'&baseline_context%in%d$baseline_context]
 contexts<-unique(o$baseline_context);ci<-match(o$baseline_context,contexts);x<-suppressWarnings(H5ADMatrix(f));lib<-DelayedMatrixStats::colSums2(x,useNames=FALSE);sumprof<-matrix(0,length(genes),length(contexts));n<-integer(length(contexts));idxg<-match(genes,feature_lists[[k]])
 for(start in seq.int(1,nrow(o),by=1000L)){
  ii<-start:min(start+999L,nrow(o));idx<-o$column_index[ii];good<-lib[idx]>0;idx<-idx[good];ii<-ii[good];b<-as(x[idxg,idx,drop=FALSE],'dgCMatrix');b@x<-log1p(b@x*rep(1e4/lib[idx],diff(b@p)));des<-sparseMatrix(i=seq_along(ii),j=ci[ii],x=1,dims=c(length(ii),length(contexts)));sumprof<-sumprof+as.matrix(b%*%des);n<-n+tabulate(ci[ii],nbins=length(contexts))
 }
 profiles[[k]]<-sweep(sumprof,2,n,'/');rownames(profiles[[k]])<-genes;colnames(profiles[[k]])<-contexts
 metadata[[k]]<-merge(data.table(baseline_context=contexts,n_baseline=n),unique(o[,.(baseline_context,cell_line,experiment,time_hours)]),by='baseline_context',sort=FALSE)
 membership[[k]]<-o[,.(source_file=basename(f),cell_id,baseline_context)]
 cat(basename(f),'contexts',length(contexts),'cells',sum(n),'\n')
}
pr<-do.call(cbind,profiles);meta<-rbindlist(metadata);stopifnot(!anyDuplicated(colnames(pr)),setequal(meta$baseline_context,colnames(pr)),setequal(d$baseline_context,colnames(pr)))
check<-merge(unique(d[,.(baseline_context,baseline_cells)]),meta,by='baseline_context');stopifnot(all(check$baseline_cells==check$n_baseline))
saveRDS(list(expression=pr,metadata=meta,modeling=d,normalization='per-cell log1p CP10k then context mean',source_files=basename(files)),file.path(out,'mixseq_raw_baseline_profiles.rds'),compress='xz')
fwrite(rbindlist(membership),file.path(out,'baseline_cell_context_membership.tsv.gz'),sep='\t');fwrite(check,file.path(out,'baseline_count_verification.tsv'),sep='\t');cat('All baseline counts match original endpoints; normalized profiles saved.\n')
