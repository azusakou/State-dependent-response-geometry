#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(args[1]);study<-args[2]
out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'));verify_locks(out)
dir.create(file.path(out,'scored'),showWarnings=FALSE);dest<-file.path(out,'scored')
started<-Sys.time();cat('Started',study,format(started,tz='UTC',usetz=TRUE),'\n');flush.console()
prep<-readRDS(file.path(out,'metadata_preparation.rds'));universe<-prep$universe;spec<-readRDS(file.path(out,'scoring_spec.rds'));real_spec<-spec[seq_len(21L)];null_names<-names(spec)[-(seq_len(21L))]
real_genes<-sort(unique(unlist(lapply(real_spec,function(s)c(s$positive,s$negative)))));study_name<-study;targets<-fread(file.path(out,'fixed_study_gene_rosters.tsv'))[study==study_name,gene]
metadata<-get_metadata(study,project,out,prep);n<-nrow(metadata);ci<-which(metadata$is_control);nc<-length(ci);ctrl_map<-integer(n);ctrl_map[ci]<-seq_len(nc)
units<-metadata[is_control==FALSE,.(n_target=.N),by=.(gene,context,gem)];setorder(units,context,gene,gem);units[,unit_index:=.I]
unit_code<-match(paste(metadata$gene,metadata$context,metadata$gem,sep='|'),paste(units$gene,units$context,units$gem,sep='|'))
unit_sums<-matrix(0,nrow(units),length(spec),dimnames=list(NULL,names(spec)))
control_base<-matrix(NA_real_,nc,length(spec),dimnames=list(NULL,names(spec)))
real_fixed<-real_exact<-matrix(NA_real_,n,21L,dimnames=list(NULL,names(real_spec)))
raw_real<-matrix(NA_real_,n,length(real_genes),dimnames=list(NULL,real_genes))
control_targets<-matrix(NA_real_,nc,length(targets),dimnames=list(NULL,targets))
for(nm in c('common_library','common_detected_genes','matrix_library','assigned_target_count'))metadata[,(nm):=NA_real_]
rawmap<-prep$rawmaps[[study]]$gene;allgenes<-unique(rawmap)
coefficient<-make_coeff(universe,spec);real_coefficient<-make_coeff(real_genes,real_spec)
source_path<-prep$paths[[study]];f<-hdf5r::H5File$new(source_path,'r')
common_raw<-which(rawmap%in%universe);aggregation<-sparseMatrix(i=match(rawmap[common_raw],universe),j=seq_along(common_raw),x=1,dims=c(length(universe),length(common_raw)))
frang_counts<-NULL
if(study=='Frangieh'){
 # CSC storage is read once by gene, retaining only selected cells and common/target columns.
 ptr<-as.numeric(f[['X']][['indptr']][]);stopifnot(length(ptr)==length(rawmap)+1L)
 gene_wanted<-unique(c(universe,intersect(targets,allgenes)));cell_map<-integer(max(metadata$source_row));cell_map[metadata$source_row]<-seq_len(n)
 ip<-jp<-xp<-vector('list',length(gene_wanted))
 for(k in seq_along(gene_wanted)){
  j<-match(gene_wanted[k],rawmap);stopifnot(!is.na(j));lo<-ptr[j]+1;hi<-ptr[j+1L]
  if(hi>=lo){rr<-as.integer(f[['X']][['indices']][lo:hi])+1L;vv<-as.numeric(f[['X']][['data']][lo:hi]);keep<-rr<=length(cell_map);ii<-cell_map[rr[keep]];vv<-vv[keep];good<-ii>0L;ip[[k]]<-ii[good];jp[[k]]<-rep.int(k,sum(good));xp[[k]]<-vv[good]}
  if(k%%1000L==0L){cat('Frangieh columns',k,'/',length(gene_wanted),'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs'))),'s\n');flush.console()}
 }
 frang_counts<-sparseMatrix(i=unlist(jp,use.names=FALSE),j=unlist(ip,use.names=FALSE),x=unlist(xp,use.names=FALSE),dims=c(length(gene_wanted),n),dimnames=list(gene_wanted,NULL));rm(ip,jp,xp);gc(FALSE)
}else if(study=='Adamson'){
 ptr<-as.numeric(f[['X']][['indptr']][]);stopifnot(length(ptr)>length(rawmap))
}
for(st in seq.int(1L,if(study=='Frangieh')n else max(metadata$source_row),by=2500L)){
 en<-min(st+2499L,if(study=='Frangieh')n else max(metadata$source_row));jj<-if(study=='Frangieh')st:en else which(metadata$source_row>=st&metadata$source_row<=en);if(!length(jj))next
 if(study=='Frangieh'){
  full<-frang_counts[,jj,drop=FALSE];b<-as.matrix(full[universe,,drop=FALSE]);lib_full<-metadata$UMI_count[jj]
  targetfun<-function(g,kk){ix<-match(g,rownames(full));if(is.na(ix))rep(NA_real_,length(kk))else as.numeric(full[ix,kk])}
 }else{
  if(study=='Adamson'){
   pp<-ptr[st:(en+1L)];ix<-seq.int(pp[1]+1,tail(pp,1));rr<-rep(seq_len(en-st+1L),diff(pp));cc<-as.integer(f[['X']][['indices']][ix])+1L;vv<-as.numeric(f[['X']][['data']][ix]);full<-t(sparseMatrix(i=rr,j=cc,x=vv,dims=c(en-st+1L,length(rawmap))))
  }else full<-f[['X']][,st:en]
  within<-metadata$source_row[jj]-st+1L;full<-full[,within,drop=FALSE]
  lib_full<-colSums(full);b<-as.matrix(aggregation%*%full[common_raw,,drop=FALSE])
  targetfun<-function(g,kk){ix<-which(rawmap==g);if(!length(ix))rep(NA_real_,length(kk))else as.numeric(colSums(full[ix,kk,drop=FALSE]))}
 }
 stopifnot(all(is.finite(b)),all(b>=0),max(abs(b-round(b)))<1e-6)
 L<-colSums(b);stopifnot(all(L>0));set(metadata,i=jj,j='common_library',value=L);set(metadata,i=jj,j='common_detected_genes',value=colSums(b>0));set(metadata,i=jj,j='matrix_library',value=lib_full)
 g<-metadata$gene[jj];target_count<-numeric(length(jj));target_count[]<-NA_real_
 for(gene in unique(g[!metadata$is_control[jj]])){kk<-which(g==gene);target_count[kk]<-targetfun(gene,kk)}
 set(metadata,i=jj,j='assigned_target_count',value=target_count)
 lg<-log1p(sweep(b,2L,L,'/')*1e4);fixed<-score_logs(lg,universe,spec,coefficient,g)
 real_fixed[jj,]<-fixed[,seq_len(21L),drop=FALSE]
 br<-b[match(real_genes,universe),,drop=FALSE];raw_real[jj,]<-t(br)
 Lminus<-L;has_target<-g%in%universe&!metadata$is_control[jj];Lminus[has_target]<-L[has_target]-target_count[has_target];stopifnot(all(is.finite(Lminus)&Lminus>0))
 lg_exact<-log1p(sweep(br,2L,Lminus,'/')*1e4);real_exact[jj,]<-score_logs(lg_exact,real_genes,real_spec,real_coefficient,g)
 cc<-which(metadata$is_control[jj]);if(length(cc)){
  cr<-ctrl_map[jj[cc]];control_base[cr,]<-fixed[cc,,drop=FALSE]
  for(gene in targets)control_targets[cr,gene]<-targetfun(gene,cc)
 }
 ti<-which(!metadata$is_control[jj]);if(length(ti)){codes<-unit_code[jj[ti]];su<-rowsum(fixed[ti,,drop=FALSE],codes,reorder=TRUE);ui<-as.integer(rownames(su));unit_sums[ui,]<-unit_sums[ui,,drop=FALSE]+su}
 if(st%%25000L==1L){cat(study,'cells through',en,'retained',sum(is.finite(metadata$common_library)),'/',n,'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s\n');flush.console()}
}
f$close_all();stopifnot(all(is.finite(real_exact)),all(is.finite(real_fixed)),all(is.finite(control_base)),all(is.finite(metadata$common_library)),sum(units$n_target)==sum(!metadata$is_control))
cache<-list(study=study,source_path=source_path,source_bytes=file.info(source_path)$size,universe=universe,real_spec=real_spec,spec=spec,targets=targets,metadata=metadata,real_scores_exact=real_exact,real_scores_fixed=real_fixed,raw_real_union=as(raw_real,'dgCMatrix'),control_target_raw=control_targets,control_cells=metadata[ci],control_rows=ci,control_fixed_base=control_base,unit_metadata=units,unit_fixed_score_sums=unit_sums,anchor_type=if(study%in%c('K562','RPE1'))'matched_GEM'else'global_within_context',membership_manifest_sha256=digest(file.path(out,'membership_lock_manifest.tsv'),file=TRUE,algo='sha256'),scoring_manifest_sha256=digest(file.path(out,'scoring_lock_manifest.tsv'),file=TRUE,algo='sha256'))
saveRDS(cache,file.path(dest,paste0(study,'_cell_cache.rds')),compress='gzip')
fwrite(metadata,file.path(dest,paste0(study,'_scored_cell_metadata.tsv.gz')),sep='\t')
fwrite(metadata[,.(n_cells=.N,technical_groups=uniqueN(gem),mean_common_library=mean(common_library),mean_common_detected_genes=mean(common_detected_genes),mean_full_library=mean(matrix_library)),by=.(study,context,gene,is_control,technical_group_type,genuine_GEM_metadata)],file.path(dest,paste0(study,'_measurement_accounting.tsv')),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(dest,paste0(study,'_extraction_sessionInfo.txt')))
cat('Completed',study,'cells',n,'controls',nc,'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'s\n')
