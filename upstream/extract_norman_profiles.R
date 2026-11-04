#!/usr/bin/env Rscript
source("code/common.R")
# Read-only aggregation of public Norman cells into condition and GEM means.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(HDF5Array);library(DelayedMatrixStats);library(hdf5r)})
source('upstream/h5ad_annotation.R')
out <- 'raw_derived';dir.create(out,recursive=TRUE,showWarnings=FALSE)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1)stop('Supply the path to the source Norman h5ad object.');input<-args[1]
features <- read_norman_input()$features
x <- H5ADMatrix(input); selected <- na.omit(match(features,rownames(x))); features <- rownames(x)[selected]
h5 <- H5File$new(input,mode='r');guides<-read_h5ad_obs(h5,'guide_identity');gem<-as.integer(read_h5ad_obs(h5,'gemgroup'));good<-tolower(as.character(read_h5ad_obs(h5,'good_coverage')));h5$close_all()
qc <- is.na(good)|good%in%c('true','1','t');cond<-sub('__.*$','',guides);levels<-unique(cond);code<-match(cond,levels)
neg<-lengths(regmatches(levels,gregexpr('NegCtrl',levels)));is_pair<-grepl('_',levels,fixed=TRUE)&neg==0;is_single<-neg==1;is_control<-neg>=2
single_gene<-vapply(strsplit(levels,'_',fixed=TRUE),function(v){v<-v[!grepl('^NegCtrl',v)];if(length(v)==1L)v else NA_character_},character(1))
total<-DelayedMatrixStats::colSums2(x,useNames=FALSE);nc<-length(levels);ng<-max(gem,na.rm=TRUE);nf<-length(features)
sums<-replicate(2,matrix(0,nf,nc),simplify=FALSE);sums_g<-replicate(2,matrix(0,nf,nc*ng),simplify=FALSE);counts<-integer(nc);counts_g<-integer(nc*ng)
for(start in seq.int(1,ncol(x),by=2000L)){
 idx<-start:min(start+1999L,ncol(x));idx<-idx[qc[idx]&is.finite(total[idx])&total[idx]>0];if(!length(idx))next
 raw<-as(x[selected,idx,drop=FALSE],'dgCMatrix');cidx<-code[idx];des<-sparseMatrix(i=seq_along(idx),j=cidx,x=1,dims=c(length(idx),nc));counts<-counts+tabulate(cidx,nbins=nc)
 vg<-is.finite(gem[idx])&gem[idx]>=1&gem[idx]<=ng;cg<-cidx[vg]+(gem[idx][vg]-1L)*nc;desg<-sparseMatrix(i=seq_len(sum(vg)),j=cg,x=1,dims=c(sum(vg),nc*ng));counts_g<-counts_g+tabulate(cg,nbins=nc*ng)
 for(k in 1:2){b<-raw;b@x<-log1p(b@x*rep(c(1e4,1e6)[k]/total[idx],diff(b@p)));sums[[k]]<-sums[[k]]+as.matrix(b%*%des);sums_g[[k]]<-sums_g[[k]]+as.matrix(b[,vg,drop=FALSE]%*%desg)}
}
inputs<-lapply(1:2,function(k){means<-sweep(sums[[k]],2,pmax(counts,1),'/');meansg<-sweep(sums_g[[k]],2,pmax(counts_g,1),'/');ctrl<-rowSums(sums[[k]][,is_control,drop=FALSE])/sum(counts[is_control]);delta<-sweep(means,1,ctrl,'-');dg<-array(NA_real_,c(nf,nc,ng),dimnames=list(features,levels,paste0('gem',1:ng)))
for(g in 1:ng){cols<-seq_len(nc)+(g-1)*nc;ctrcols<-which(is_control)+(g-1)*nc;ctrlg<-rowSums(sums_g[[k]][,ctrcols,drop=FALSE])/sum(counts_g[ctrcols]);dg[,,g]<-sweep(meansg[,cols,drop=FALSE],1,ctrlg,'-');dg[,counts_g[cols]==0,g]<-NA_real_}
list(delta=delta,delta_gem=dg,normalization=c('log1p_CP10k','log1p_CPM')[k])})
common<-list(features=features,conditions=levels,n_cells=counts,n_condition_gem=matrix(counts_g,nc,ng),is_pair=is_pair,is_single=is_single,is_control=is_control,single_gene=single_gene)
primary<-c(common,inputs[[1]])
saveRDS(primary,file.path(out,'norman_condition_gem_input.rds'),compress='gzip')
saveRDS(c(common,inputs[[2]]),file.path(out,'norman_CPM_sensitivity_input.rds'),compress='gzip')
fwrite(data.table(gene=features),file.path(out,'norman_response_genes.tsv'),sep='\t')
fwrite(data.table(cell_id=colnames(x),condition_id=cond,gem_group=gem,good_coverage=qc,positive_library=is.finite(total)&total>0),file.path(out,'norman_cell_assignments.tsv.gz'),sep='\t')
cat('Extracted',nf,'genes,',nc,'conditions,',ng,'technical GEM groups.\n')
