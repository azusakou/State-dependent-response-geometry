#!/usr/bin/env Rscript
# Reconstruct identically defined per-cell scores from read-only count matrices.
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(args[1]); out <- args[2]
.libPaths(c(file.path(project,'.Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(yaml);library(Matrix);library(HDF5Array);library(hdf5r)})
dir.create(out, recursive=TRUE, showWarnings=FALSE)
annot <- function(f,g,n) {
  node <- f[[g]][[n]]
  if(inherits(node,'H5Group')) {
    codes <- as.integer(node[['codes']][]); cats <- node[['categories']][]
    ans <- rep(NA_character_,length(codes));ok<-codes>=0;ans[ok]<-cats[codes[ok]+1L];return(ans)
  }
  ans <- node[]
  if('__categories'%in%names(f[[g]]) && n%in%names(f[[g]][['__categories']])) {
    cats <- f[[g]][['__categories']][[n]][];ans<-cats[as.integer(ans)+1L]
  }
  ans
}
repl <- '/Volumes/wang/paper6yhn/raw/replogle_k562_essential/K562_essential_raw_singlecell_01.h5ad'
norman <- '/Volumes/wang/paper4rfv/data/processed/h5ad/core/norman_2019_combinatorial.h5ad'
mixdir <- '/Volumes/wang/paper6yhn/processed/mixseq'
mixfiles<-list.files(mixdir,pattern='\\.qc\\.h5ad$',full.names=TRUE)
mixfiles<-mixfiles[!grepl('/\\._',mixfiles)]
spec<-yaml::read_yaml(file.path(project,'paper19opp-revision1/GitHub/upstream/program_gene_sets.yaml'))$axes
hallmark_file<-'/Users/wang/Documents/paper6yhn/resource/networks/hallmark_human.tsv.gz'
hallmark<-fread(hallmark_file)
pathways<-c('UNFOLDED_PROTEIN_RESPONSE','E2F_TARGETS','G2M_CHECKPOINT','P53_PATHWAY','TNFA_SIGNALING_VIA_NFKB','INTERFERON_GAMMA_RESPONSE','MTORC1_SIGNALING','OXIDATIVE_PHOSPHORYLATION','HYPOXIA')
stopifnot(all(pathways%in%hallmark$source))
for(a in pathways)spec[[paste0('hallmark_',a)]]<-list(positive_genes=unique(hallmark[source==a,target]),negative_genes=character())
files<-c(repl,norman,mixfiles)
feature_lists<-lapply(files,function(p){f<-hdf5r::H5File$new(p,'r');on.exit(f$close_all());
 if('gene_name'%in%names(f[['var']]))as.character(annot(f,'var','gene_name'))else if('gene_symbol'%in%names(f[['var']]))as.character(annot(f,'var','gene_symbol'))else if('gene'%in%names(f[['var']]))as.character(annot(f,'var','gene'))else as.character(annot(f,'var','_index'))})
external_paths<-file.path(project,c('paper19opp-revision2/reconstruction/candidate_validation/adamson_measured_gene_symbols.txt','paper19opp-revision2/reconstruction/external_resources/frangieh_measured_gene_symbols.txt'))
stopifnot(all(file.exists(external_paths)))
common<-sort(Reduce(intersect,c(feature_lists,lapply(external_paths,readLines))))
stopifnot(length(common)>1000,!anyDuplicated(common))
fwrite(data.table(gene=common),file.path(out,'final_gene_universe.tsv'),sep='\t')
coverage<-rbindlist(lapply(names(spec),function(a)rbindlist(lapply(c('positive_genes','negative_genes'),function(direction){g<-spec[[a]][[direction]];if(!length(g))return(NULL);data.table(program=a,direction=direction,gene=g,in_common=g%in%common)}))))
fwrite(coverage,file.path(out,'final_module_gene_coverage.tsv'),sep='\t')
fwrite(coverage[in_common==TRUE],file.path(out,'final_axis_membership.tsv'),sep='\t')
fwrite(data.table(source_file=files,bytes=file.info(files)$size,modified=as.character(file.info(files)$mtime),n_features=lengths(feature_lists),n_duplicate_symbols=vapply(feature_lists,function(g)sum(duplicated(g)),integer(1))),file.path(out,'final_source_inventory.tsv'),sep='\t')
writeLines(c('metadata-only feature lock; no validation effects inspected',paste('Hallmark source:',hallmark_file),paste('Measured external symbols:',external_paths),paste('Frozen axes:',paste(names(spec),collapse=', '))),file.path(out,'feature_lock.txt'))
cat('COMMON_GENES',length(common),'\n');print(coverage[,.(present=sum(in_common),specified=.N),by=.(program,direction)])
rq<-as.data.table(readRDS(file.path(project,'paper19opp-revision1/GitHub/inputs/replogle_scores.rds')))
extract<-function(k){
 p<-files[k];f<-hdf5r::H5File$new(p,'r')
 if(k==1L){
  o<-rq[,.(cell_id,gene=as.character(gene),gem_group,passes_qc_final,UMI_count,n_genes_by_counts,mitopercent)]
  o[,`:=`(dataset='Replogle',context=paste0('GEM',gem_group),biological_context='K562',perturbation=gene,is_control=gene=='non-targeting',mode='CRISPRi')]
 }else if(k==2L){
  cond<-sub('__.*$','',as.character(annot(f,'obs','guide_identity')))
  neg<-lengths(regmatches(cond,gregexpr('NegCtrl',cond)))
  o<-data.table(cell_id=as.character(annot(f,'obs','_index')),perturbation=cond,gem_group=as.integer(annot(f,'obs','gemgroup')),is_control=neg>=2L,passes_qc_final=tolower(as.character(annot(f,'obs','good_coverage')))%in%c('true','1','t'))
  o[,`:=`(dataset='Norman',context=paste0('GEM',gem_group),biological_context='K562',mode='CRISPRa')]
 }else{
  o<-data.table(cell_id=as.character(annot(f,'obs','cell_id')),cell_line=as.character(annot(f,'obs','cell_line')),experiment=as.integer(annot(f,'obs','experiment')),time_hours=as.numeric(annot(f,'obs','time_hours')),task=as.character(annot(f,'obs','task')),perturbation=as.character(annot(f,'obs','perturbation')),passes_qc_final=as.logical(annot(f,'obs','passes_qc_final')))
  o[,`:=`(dataset='MIX-Seq',context=paste(experiment,ifelse(is.na(time_hours),'NA',time_hours),cell_line,task,sep='|'),biological_context=cell_line,is_control=tolower(perturbation)%in%c('dmso','untreated','sglacz','sgor2j2'),mode=ifelse(task=='chemical','drug','CRISPR'))]
 }
 f$close_all();o[,column_index:=.I];o<-o[passes_qc_final==TRUE]
 # DMSO and sgLACZ take priority when additional contemporaneous controls exist.
 o[,source_file:=basename(p)]
 x<-suppressWarnings(H5ADMatrix(p));stopifnot(nrow(x)==length(feature_lists[[k]]))
 gi<-which(feature_lists[[k]]%in%common)
 # Aggregate all columns/rows sharing a measured gene symbol before log1p.
 aggregation<-sparseMatrix(i=match(feature_lists[[k]][gi],common),j=seq_along(gi),x=1,dims=c(length(common),length(gi)))
 idx<-lapply(names(spec),function(a){list(pos=which(common%in%spec[[a]]$positive_genes),neg=which(common%in%spec[[a]]$negative_genes))});names(idx)<-names(spec)
 for(a in names(spec))o[,(paste0('score_',a)):=NA_real_]
 o[,`:=`(common_library=NA_real_,common_detected_genes=NA_integer_)]
 control_sum<-numeric(length(common));control_n<-0L
 for(st in seq.int(1,nrow(o),by=1500L)){
  ii<-st:min(st+1499L,nrow(o));b<-as.matrix(aggregation%*%as.matrix(x[gi,o$column_index[ii],drop=FALSE]));lib<-colSums(b)
  o$common_detected_genes[ii]<-colSums(b>0)
  b<-log1p(sweep(b,2,pmax(lib,1),'/')*1e4)
  o$common_library[ii]<-lib
  co<-which(o$is_control[ii]&lib>0);if(length(co)){control_sum<-control_sum+rowSums(b[,co,drop=FALSE]);control_n<-control_n+length(co)}
  for(a in names(spec)){
   pos<-colMeans(b[idx[[a]]$pos,,drop=FALSE]);neg<-if(length(idx[[a]]$neg))colMeans(b[idx[[a]]$neg,,drop=FALSE])else 0
   set(o,i=ii,j=paste0('score_',a),value=pos-neg)
  }
 }
 o<-o[is.finite(common_library)&common_library>0]
 o[,c('column_index','passes_qc_final'):=NULL]
 fwrite(data.table(gene=common,mean_control_log1p_common_CP10k=control_sum/max(control_n,1),n_control_cells=control_n),file.path(out,paste0('control_gene_expression_source',k,'.tsv')),sep='\t')
 cat(k,basename(p),'scored',nrow(o),'\n');flush.console();o
}
parts<-vector('list',length(files))
for(k in seq_along(files)){
 cache<-file.path(out,paste0('harmonized_score_part_',k,'.rds'))
 if(file.exists(cache))parts[[k]]<-readRDS(cache)else{parts[[k]]<-extract(k);saveRDS(parts[[k]],cache,compress='gzip')}
}
d<-rbindlist(parts,fill=TRUE)
# Match the original source-specific choice of controls, without endpoint filtering.
d[dataset=='MIX-Seq',control_priority:=ifelse(tolower(perturbation)=='dmso',1L,ifelse(tolower(perturbation)=='sglacz',1L,2L))]
d<-d[!(dataset=='MIX-Seq'&is_control&control_priority>1L)]
fwrite(d,file.path(out,'common_cell_scores.tsv.gz'),sep='\t')
saveRDS(d,file.path(out,'common_cell_scores.rds'),compress='gzip')
cat('FINAL_CELLS',nrow(d),'\n');print(d[,.(cells=.N,contexts=uniqueN(context),interventions=uniqueN(perturbation[!is_control]),control_cells=sum(is_control)),by=dataset])
writeLines(capture.output(sessionInfo()),file.path(out,'score_session_info.txt'))
