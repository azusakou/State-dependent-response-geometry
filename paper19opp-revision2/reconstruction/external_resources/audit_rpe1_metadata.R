#!/usr/bin/env Rscript
# Read identifiers and QC metadata only. Never access RPE1 expression matrix values.
suppressPackageStartupMessages({library(hdf5r);library(data.table)})
out<-'paper19opp-revision2/reconstruction/external_resources';p<-'/Volumes/wang/paper6yhn/raw/replogle_rpe1/rpe1_raw_singlecell_01.h5ad'
f<-H5File$new(p,'r')
readmeta<-function(group,name){g<-f[[group]];v<-g[[name]][];if('__categories'%in%g$ls()$name && name%in%g[['__categories']]$ls()$name){lev<-g[['__categories']][[name]][];ans<-rep(NA_character_,length(v));ok<-v>=0&v<length(lev);ans[ok]<-lev[v[ok]+1L];v<-ans};v}
genes<-readmeta('var','gene_name');ids<-readmeta('var','gene_id')
universe<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_gene_universe.tsv')$gene
axes<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_axis_membership.tsv');upr<-axes[program=='hallmark_UNFOLDED_PROTEIN_RESPONSE'&in_common==TRUE,gene]
md<-data.table(source_row=seq_along(readmeta('obs','cell_barcode')),cell_id=readmeta('obs','cell_barcode'),gene=readmeta('obs','gene'),gem_group=readmeta('obs','gem_group'),sgID_AB=readmeta('obs','sgID_AB'),UMI_count=readmeta('obs','UMI_count'),mitopercent=readmeta('obs','mitopercent'))
md[,is_control:=gene=='non-targeting']
fwrite(md,file.path(out,'rpe1_metadata_only_cells.tsv.gz'),sep='\t')
fwrite(md[,.(cells=.N,gem_groups=uniqueN(gem_group),guide_constructs=uniqueN(sgID_AB)),by=.(gene,is_control)],file.path(out,'rpe1_metadata_only_gene_counts.tsv'),sep='\t')
fwrite(md[,.(cells=.N,control_cells=sum(is_control)),by=gem_group],file.path(out,'rpe1_metadata_only_gem_counts.tsv'),sep='\t')
fwrite(data.table(feature_row=seq_along(genes),gene=genes,gene_id=ids),file.path(out,'rpe1_metadata_only_feature_mapping.tsv'),sep='\t')
fwrite(data.table(gene=universe,measured=universe%in%genes,is_upr=universe%in%upr),file.path(out,'rpe1_metadata_only_universe_coverage.tsv'),sep='\t')
obs_names<-f[['obs']]$ls()$name
writeLines(c(paste('source',p),paste('source_bytes',file.info(p)$size),paste('source_cells',nrow(md)),paste('measured_feature_rows',length(genes)),paste('unique_symbols',uniqueN(genes)),paste('duplicate_symbol_rows',sum(duplicated(genes))),paste('common_universe_present',sum(universe%in%genes),'of',length(universe)),paste('UPR_members_present',sum(upr%in%genes),'of',length(upr)),paste('controls',sum(md$is_control)),paste('target_genes',uniqueN(md[is_control==FALSE]$gene)),paste('GEM_groups',uniqueN(md$gem_group)),paste('QC_flag_columns',paste(grep('qc|pass|filter|doublet',obs_names,ignore.case=TRUE,value=TRUE),collapse=',')),paste('observed_metadata_columns',paste(obs_names,collapse=',')),'RPE1 expression matrix values and var mean/std/cv/fano values were not accessed. QC flags were not invented; UMI_count and mitopercent are available, while gene/construct/GEM counts are metadata-only eligibility information.'),file.path(out,'rpe1_metadata_only_summary.txt'))
f$close_all();cat(readLines(file.path(out,'rpe1_metadata_only_summary.txt')),sep='\n')
