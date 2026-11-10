#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(hdf5r);library(Matrix);library(yaml);library(digest)})
a<-commandArgs(TRUE);inp<-if(length(a))a[1]else'/Volumes/wang/paper4rfv/data/processed/h5ad/core_unpacked/adamson_2016_upr__GSM2406681_10X010.h5ad';out<-if(length(a)>1)a[2]else'paper19opp-revision2/reconstruction/candidate_validation';dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(out,'discovery_lock_hashes.tsv')))
lock<-fread(file.path(out,'discovery_lock_hashes.tsv'));stopifnot(all(vapply(file.path(out,lock$file),digest,character(1),file=TRUE,algo='sha256')==lock$sha256))
membership<-fread(file.path(out,'final_axis_membership.tsv'));sets<-lapply(unique(membership$program),function(ax){list(positive_genes=membership[program==ax&direction=='positive_genes'&in_common==TRUE,gene],negative_genes=membership[program==ax&direction=='negative_genes'&in_common==TRUE,gene])});names(sets)<-unique(membership$program)
h<-H5File$new(inp,'r');ann<-function(gr,f){nd<-h[[gr]][[f]];if(inherits(nd,'H5Group')){co<-nd[['codes']][];cc<-nd[['categories']][];z<-rep(NA_character_,length(co));z[co>=0]<-cc[co[co>=0]+1];z}else nd[]}
genes<-as.character(ann('var','_index'));meta<-data.table(cell_id=ann('obs','_index'),gene=ann('obs','perturbation'),guide=ann('obs','guide_identity'));meta[,context:='Adamson_K562_UPR'];meta[,is_control:=gene=='63(mod)'];meta[,is_alternative_control:=gene=='Gal4-4(mod)']
universe<-fread(file.path(out,'final_gene_universe.tsv'))$gene;commonix<-which(genes%in%universe);stopifnot(length(commonix)==7848,length(unique(genes[commonix]))==7848);ng<-length(genes);n<-nrow(meta);ptr<-as.numeric(h[['X']][['indptr']][])
# Scores use precisely the original positive-minus-negative log1pCP10k rule.
coverage<-rbindlist(lapply(names(sets),function(axis)rbindlist(lapply(c('positive_genes','negative_genes'),function(sign){wanted<-unlist(sets[[axis]][[sign]]);if(!length(wanted))return(NULL);data.table(axis=rep(axis,length(wanted)),direction=rep(sign,length(wanted)),gene=wanted,present=wanted%in%genes)}))))
fwrite(coverage,file.path(out,'adamson_program_feature_coverage_harmonized.tsv'),sep='\t')
targetgenes<-sort(unique(meta$gene[meta$gene%in%fread(file.path(out,'locked_discovery_ranking.tsv'))$gene]));targetix<-match(targetgenes,genes)
control_sum<-control_nz<-numeric(length(targetgenes));control_n<-0L
outrows<-list()
for(st in seq(1,n,by=2000L)){
 en<-min(n,st+1999L);pp<-ptr[st:(en+1)];ind<-seq.int(pp[1]+1,tail(pp,1));jj<-as.integer(h[['X']][['indices']][ind])+1L;vv<-as.numeric(h[['X']][['data']][ind]);ii<-rep(seq_len(en-st+1L),diff(pp))
 mat<-sparseMatrix(i=ii,j=jj,x=vv,dims=c(en-st+1L,ng));lib_full<-Matrix::rowSums(mat);lib<-Matrix::rowSums(mat[,commonix,drop=FALSE]);nnz<-Matrix::rowSums(mat>0)
 norm<-mat;norm@x<-log1p(norm@x/pmax(lib[norm@i+1L],1)*1e4)
 z<-copy(meta[st:en]);z[,`:=`(library_size=lib,full_library_size=lib_full,detected_genes=nnz)]
 for(axis in names(sets)){
  p<-match(unlist(sets[[axis]]$positive_genes),genes);p<-p[!is.na(p)];q<-match(unlist(sets[[axis]]$negative_genes),genes);q<-q[!is.na(q)]
  z[,(paste0('score_',axis)):=if(length(p))Matrix::rowMeans(norm[,p,drop=FALSE])else NA_real_]
  if(length(q))z[,(paste0('score_',axis)):=get(paste0('score_',axis))-Matrix::rowMeans(norm[,q,drop=FALSE])]
  allnorm<-mat;allnorm@x<-log1p(allnorm@x/pmax(lib_full[allnorm@i+1L],1)*1e4);sv<-if(length(p))Matrix::rowMeans(allnorm[,p,drop=FALSE])else rep(NA_real_,nrow(z));if(length(q))sv<-sv-Matrix::rowMeans(allnorm[,q,drop=FALSE]);z[,(paste0('full_library_score_',axis)):=sv]
 }
 self<-match(z$gene,genes);good<-which(!is.na(self));z[,`:=`(target_raw_count=NA_real_,target_logcp10k=NA_real_)];z[good,target_raw_count:=as.numeric(mat[cbind(good,self[good])])];z[good,target_logcp10k:=as.numeric(norm[cbind(good,self[good])])]
 ctl<-which(z$is_control);if(length(ctl)){present<-which(!is.na(targetix));control_sum[present]<-control_sum[present]+Matrix::colSums(norm[ctl,targetix[present],drop=FALSE]);control_nz[present]<-control_nz[present]+Matrix::colSums(mat[ctl,targetix[present],drop=FALSE]>0);control_n<-control_n+length(ctl)}
 outrows[[length(outrows)+1L]]<-z
 cat('Scored',en,'of',n,'\n')
}
h$close_all();z<-rbindlist(outrows);fwrite(z,file.path(out,'adamson_cell_scores_harmonized.tsv.gz'),sep='\t')
fwrite(data.table(gene=targetgenes,target_measured=!is.na(targetix),control_mean_logcp10k=control_sum/max(control_n,1),control_detected_fraction=control_nz/max(control_n,1),control_cells=control_n),file.path(out,'adamson_control_target_expression_harmonized.tsv'),sep='\t')
fwrite(z[,.(cells=.N,guide_constructs=uniqueN(guide)),by=.(gene,context,is_control,is_alternative_control)],file.path(out,'adamson_scored_counts_harmonized.tsv'),sep='\t')
fwrite(data.table(source=normalizePath(inp),bytes=file.info(inp)$size,sha256=digest(inp,file=TRUE,algo='sha256'),accession='GSE90546/GSM2406681',source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681'),file.path(out,'adamson_source_provenance_harmonized.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'extraction_harmonized_sessionInfo.txt'))
