#!/usr/bin/env Rscript
project<-normalizePath(commandArgs(TRUE)[1]);out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'));verify_locks(out);dest<-file.path(out,'canonical_Hallmark94')
lock<-fread(file.path(dest,'technical_addendum_lock.tsv'));stopifnot(all(vapply(file.path(dest,lock$file),digest,character(1),file=TRUE,algo='sha256')==lock$sha256))
aliases<-c('H2AFX','IARS','TARS','TTC37');prep<-readRDS(file.path(out,'metadata_preparation.rds'));account<-list()
for(study_name in c('K562','RPE1','Adamson','Frangieh')){
 started<-Sys.time();cat('Narrow four-gene extraction',study_name,format(started,tz='UTC',usetz=TRUE),'\n');flush.console()
 cache<-readRDS(file.path(out,'scored',paste0(study_name,'_cell_cache.rds')));md<-cache$metadata;genes<-prep$rawmaps[[study_name]]$gene;stopifnot(all(aliases%in%genes),!any(aliases%in%colnames(cache$raw_real_union)))
 source_path<-prep$paths[[study_name]];f<-hdf5r::H5File$new(source_path,'r');additional<-matrix(0,nrow(md),4L,dimnames=list(md$cell_id,aliases))
 if(study_name%in%c('K562','RPE1')){
  ix<-which(genes%in%aliases);allraw<-f[['X']][ix,];stopifnot(nrow(allraw)==length(ix),ncol(allraw)>=max(md$source_row))
  for(g in aliases)additional[,g]<-colSums(allraw[genes[ix]==g,md$source_row,drop=FALSE])
  native_cells<-ncol(allraw)
 }else if(study_name=='Adamson'){
  ptr<-as.numeric(f[['X']][['indptr']][]);native_cells<-length(ptr)-1L;source_map<-integer(native_cells);source_map[md$source_row]<-seq_len(nrow(md));gi<-which(genes%in%aliases)
  for(st in seq.int(1L,native_cells,by=2500L)){
   en<-min(st+2499L,native_cells);pp<-ptr[st:(en+1L)];ix<-seq.int(pp[1]+1,tail(pp,1));cc<-as.integer(f[['X']][['indices']][ix])+1L;vv<-as.numeric(f[['X']][['data']][ix]);rr<-rep.int(seq.int(st,en),diff(pp));keep<-cc%in%gi&source_map[rr]>0L
   if(any(keep)){ii<-source_map[rr[keep]];jj<-match(genes[cc[keep]],aliases);additional[cbind(ii,jj)]<-vv[keep]}
  }
 }else{
  ptr<-as.numeric(f[['X']][['indptr']][]);native_cells<-218331L
  for(g in aliases){
   j<-match(g,genes);lo<-ptr[j]+1;hi<-ptr[j+1L];v<-numeric(native_cells)
   if(hi>=lo)v[as.integer(f[['X']][['indices']][lo:hi])+1L]<-as.numeric(f[['X']][['data']][lo:hi])
   additional[,g]<-v[md$source_row]
  }
 }
 f$close_all();stopifnot(all(is.finite(additional)),all(additional>=0),max(abs(additional-round(additional)))<1e-6)
 shared<-intersect(aliases,colnames(cache$control_target_raw));err<-if(length(shared))max(abs(additional[cache$control_rows,shared,drop=FALSE]-cache$control_target_raw[,shared,drop=FALSE]))else NA_real_;stopifnot(!is.finite(err)||err==0)
 result<-list(study=study_name,source_path=source_path,source_bytes=file.info(source_path)$size,cell_id=md$cell_id,source_row=md$source_row,raw_counts=additional,shared_control_targets=shared,maximum_shared_control_count_error=err,addendum_manifest_sha256=digest(file.path(dest,'technical_addendum_lock.tsv'),file=TRUE,algo='sha256'))
 saveRDS(result,file.path(dest,paste0(study_name,'_additional_raw_counts.rds')),compress='gzip')
 account[[study_name]]<-data.table(study=study_name,gene=aliases,raw_feature_rows=vapply(aliases,function(g)paste(which(genes==g),collapse=';'),character(1)),native_cells=native_cells,retained_cells=nrow(md),shared_control_target_check=aliases%in%shared,max_shared_control_error=err,source_path=source_path,source_bytes=file.info(source_path)$size,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
 cat('Complete',study_name,'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),2),'s; overlapQAerror',err,'\n');flush.console();rm(cache,additional,result);gc(FALSE)
}
fwrite(rbindlist(account),file.path(dest,'narrow_extraction_accounting.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(dest,'narrow_extraction_sessionInfo.txt'))
