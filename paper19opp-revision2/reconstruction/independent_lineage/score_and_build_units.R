#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix)})
out<-'paper19opp-revision2/reconstruction/independent_lineage';started<-Sys.time()
cat('Scoring started',format(started,tz='UTC'),'UTC\n');flush.console()
stopifnot(file.exists(file.path(out,'locked_analysis_protocol.sha256')))
counts<-readRDS(file.path(out,'primary_counts.rds'));md<-fread(file.path(out,'primary_cell_metadata.tsv.gz'))
stopifnot(identical(colnames(counts),md$cell_id),!anyDuplicated(rownames(counts)),all(counts@x>=0),all(is.finite(counts@x)),max(abs(counts@x-round(counts@x)))<1e-6)
full_lib<-Matrix::colSums(counts);stopifnot(max(abs(full_lib-md$nCount_RNA))<1e-5)
universe<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_gene_universe.tsv')$gene
present<-intersect(universe,rownames(counts));stopifnot(length(present)==7848L)
writeLines(present,file.path(out,'used_gene_universe.txt'))
idx<-match(present,rownames(counts));z<-counts[idx,,drop=FALSE];lib<-Matrix::colSums(z);stopifnot(all(lib>0))
z<-z%*%Matrix::Diagonal(x=1e4/lib);z@x<-log1p(z@x)
axes<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_axis_membership.tsv')[in_common==TRUE]
coverage<-list();features<-character()
for(a in unique(axes$program)){
 pos<-axes[program==a&direction=='positive_genes',gene];neg<-axes[program==a&direction=='negative_genes',gene]
 stopifnot(all(c(pos,neg)%in%present),length(pos)>0)
 v<-Matrix::colMeans(z[match(pos,present),,drop=FALSE]);if(length(neg))v<-v-Matrix::colMeans(z[match(neg,present),,drop=FALSE])
 nm<-paste0('score_',a);md[,(nm):=v];features<-c(features,nm)
 coverage[[a]]<-data.table(axis=a,positive_genes=length(pos),negative_genes=length(neg))
}
rm(z,counts);gc()
md[,`:=`(lineage=Lineage,first=sub('to.*$','',OG_condition),second=ifelse(grepl('to',OG_condition),sub('^.*?to','',OG_condition),''),full_library=full_lib,common_library=lib)]
stopifnot(all(md$first%in%c('cis','cocl2','dabtram')),all(md$second%in%c('','cis','cocl2','dabtram')))
md[,valid_lineage:=grepl('^Lin[0-9]+$',lineage)]
fwrite(md[,.(cells=.N,valid_lineage_cells=sum(valid_lineage),valid_lineages=uniqueN(lineage[valid_lineage])),by=.(first,second)],file.path(out,'source_sample_accounting.tsv'),sep='\t')
valid<-md[valid_lineage==TRUE]
centroids<-valid[,c(list(n_cells=.N,mean_log_library=mean(log1p(full_library)),mean_detected_genes=mean(nFeature_RNA)),lapply(.SD,mean)),by=.(lineage,first,second),.SDcols=features]
fwrite(centroids,file.path(out,'reconstructed_lineage_states.tsv.gz'),sep='\t')
parent<-centroids[second==''];parent[,second:=NULL];setnames(parent,setdiff(names(parent),c('lineage','first')),paste0('parent_',setdiff(names(parent),c('lineage','first'))))
desc<-centroids[second!=''];setnames(desc,setdiff(names(desc),c('lineage','first','second')),paste0('desc_',setdiff(names(desc),c('lineage','first','second'))))
units<-merge(parent,desc,by=c('lineage','first'),allow.cartesian=TRUE)
stopifnot(!anyDuplicated(units[,.(lineage,first,second)]))
units[,sequence:=paste(first,second,sep='_to_')]
for(a in c('apoptosis','proliferation','stress'))units[,(paste0('shift_',a)):=get(paste0('desc_score_',a))-get(paste0('parent_score_',a))]
fwrite(units,file.path(out,'reconstructed_paired_units.tsv.gz'),sep='\t')
all_info<-rbindlist(lapply(list(primary=c(3,2),strict=c(10,5)),function(g){u<-units[parent_n_cells>=g[1]&desc_n_cells>=g[2]];u[,.(tasks=.N,lineages=uniqueN(lineage)),by=.(first,second)]}),idcol='support')
fwrite(all_info,file.path(out,'eligible_unit_accounting.tsv'),sep='\t')
fwrite(rbindlist(coverage),file.path(out,'axis_gene_coverage.tsv'),sep='\t')
fwrite(valid[,c('cell_id','lineage','first','second','full_library','common_library','nFeature_RNA',features),with=FALSE],file.path(out,'scored_lineage_cells.tsv.gz'),sep='\t')
saveRDS(list(units=units,features=features),file.path(out,'analysis_units.rds'))
writeLines(c(paste('source_cells',nrow(md)),paste('valid_lineage_cells',nrow(valid)),paste('unique_assigned_lineages',uniqueN(valid$lineage)),paste('paired_tasks',nrow(units)),paste('paired_unique_lineages',uniqueN(units$lineage)),paste('elapsed_seconds',as.numeric(difftime(Sys.time(),started,units='secs')))),file.path(out,'scoring_summary.txt'))
cat('Scoring complete in',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n');print(all_info)
