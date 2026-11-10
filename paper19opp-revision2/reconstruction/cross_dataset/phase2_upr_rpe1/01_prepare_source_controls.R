#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);project<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(project,'.Rlib'),.libPaths()));source(file.path(out,'phase2_helpers.R'))
dir.create(out,showWarnings=FALSE,recursive=TRUE)
cross<-dirname(out);external<-file.path(project,'paper19opp-revision2/reconstruction/external_resources');candidates<-file.path(project,'paper19opp-revision2/reconstruction/candidate_validation')
original_universe<-fread(file.path(cross,'final_gene_universe.tsv'))$gene
coverage<-fread(file.path(external,'rpe1_metadata_only_universe_coverage.tsv'))
universe<-sort(coverage[measured==TRUE,gene]);stopifnot(length(universe)==6780L,all(universe%in%original_universe))
original_upr<-fread(file.path(cross,'final_axis_membership.tsv'))[program=='hallmark_UNFOLDED_PROTEIN_RESPONSE',gene]
upr<-sort(intersect(original_upr,universe));stopifnot(length(original_upr)==91L,length(upr)==90L,setequal(setdiff(original_upr,upr),'NFYA'))
fwrite(data.table(gene=universe),file.path(out,'phase2_gene_universe.tsv'),sep='\t')
fwrite(data.table(gene=upr),file.path(out,'phase2_UPR_members.tsv'),sep='\t')
rq<-as.data.table(readRDS(file.path(project,'paper19opp-revision1/GitHub/inputs/replogle_scores.rds')))
rq[,source_row:=.I];rq[,gene:=as.character(gene)];rq<-rq[passes_qc_final==TRUE]
meta<-rq[,.(source_row,cell_id,gene,gem_group,UMI_count,mitopercent,n_genes_by_counts,is_control=gene=='non-targeting')]
source_gem<-meta[is_control==FALSE,.N,by=.(gene,gem_group)]
source_count<-source_gem[,.(source_cells=sum(N),source_gem_labels=.N,source_GEM_ge5=sum(N>=5)),by=gene]
rpe<-fread(file.path(external,'rpe1_metadata_only_cells.tsv.gz'));rpe<-rpe[is.finite(UMI_count)&UMI_count>0]
rpe_count<-rpe[is_control==FALSE,.(rpe_cells=.N,rpe_gem_labels=uniqueN(gem_group)),by=gene]
rpe_gem<-rpe[is_control==FALSE,.N,by=.(gene,gem_group)]
rpe_strict<-rpe_gem[,.(rpe_GEM_ge5=sum(N>=5)),by=gene]
panel<-merge(merge(source_count,rpe_count,by='gene'),rpe_strict,by='gene')
panel[,eligible_primary:=source_cells>=100&source_GEM_ge5>=3&rpe_cells>=100&rpe_gem_labels>=3]
panel[,eligible_weak_source:=source_cells>=100&source_gem_labels>=3&rpe_cells>=100&rpe_gem_labels>=3]
panel[,strict_rpe_support:=eligible_primary&rpe_GEM_ge5>=3]
fwrite(panel,file.path(out,'metadata_eligibility_accounting.tsv'),sep='\t')
eligible<-panel[eligible_primary==TRUE][order(gene)]
stopifnot(nrow(eligible)==438L,sum(panel$eligible_weak_source)==454L,sum(panel$strict_rpe_support)==390L)
fwrite(eligible,file.path(out,'phase2_eligible_genes.tsv'),sep='\t')
fwrite(panel[eligible_weak_source==TRUE&eligible_primary==FALSE],file.path(out,'weak_to_primary_exclusions_metadata_only.tsv'),sep='\t')
fwrite(meta,file.path(out,'source_QC_cells.tsv.gz'),sep='\t')
fwrite(rpe,file.path(out,'RPE1_author_roster_metadata_only.tsv.gz'),sep='\t')
fixed<-sort(unique(c(fread(file.path(candidates,'harmonized_adamson_candidates.tsv'))$gene,fread(file.path(candidates,'harmonized_frangieh_candidates.tsv'))$gene)))
target_genes<-sort(unique(c(eligible$gene,fixed)))
fwrite(data.table(gene=target_genes,in_original_fixed_external_union=target_genes%in%fixed),file.path(out,'cached_target_gene_union.tsv'),sep='\t')
source_path<-'/Volumes/wang/paper6yhn/raw/replogle_k562_essential/K562_essential_raw_singlecell_01.h5ad'
f<-hdf5r::H5File$new(source_path,'r');raw_gene<-as.character(read_annot(f,'var','gene_name'));raw_cell_gene<-as.character(read_annot(f,'obs','gene'));raw_gem<-as.integer(read_annot(f,'obs','gem_group'));f$close_all()
stopifnot(identical(raw_cell_gene[meta$source_row],meta$gene),identical(raw_gem[meta$source_row],as.integer(meta$gem_group)))
allgenes<-unique(raw_gene);aggregation<-sparseMatrix(i=match(raw_gene,allgenes),j=seq_along(raw_gene),x=1,dims=c(length(allgenes),length(raw_gene)))
x<-suppressWarnings(H5ADMatrix(source_path));stopifnot(all(universe%in%allgenes))
controls<-meta[is_control==TRUE];sumlog<-numeric(length(allgenes));detected<-sumlog;n<-0L
for(st in seq.int(1,nrow(controls),by=1200L)){
 ii<-st:min(st+1199L,nrow(controls));b<-as.matrix(aggregation%*%as.matrix(x[,controls$source_row[ii],drop=FALSE]));lib<-colSums(b[match(universe,allgenes),,drop=FALSE]);stopifnot(all(lib>0))
 sumlog<-sumlog+rowSums(log1p(sweep(b,2,lib,'/')*1e4));detected<-detected+rowSums(b>0);n<-n+length(ii)
 cat('Source NTC preparation',n,'/',nrow(controls),'\n');flush.console()
}
stats<-data.table(gene=allgenes,source_NTC_mean_logCP10k=sumlog/n,source_NTC_detection=detected/n,source_NTC_cells=n,measured_phase2=allgenes%in%universe)
fwrite(stats,file.path(out,'source_NTC_gene_statistics.tsv'),sep='\t')
pool<-stats[measured_phase2==TRUE];pool[,expression_bin:=pmin(5L,as.integer((frank(source_NTC_mean_logCP10k,ties.method='average')-1)/.N*5)+1L)];pool[,detection_bin:=pmin(5L,as.integer((frank(source_NTC_detection,ties.method='average')-1)/.N*5)+1L)]
fwrite(pool,file.path(out,'source_NTC_matching_bins.tsv'),sep='\t')
nulls<-list()
for(m in 1:100){
 set.seed(20261106L+m);used<-character()
 for(g in upr){
  ref<-pool[gene==g];available<-pool[!gene%in%c(upr,used)];available[,distance:=abs(expression_bin-ref$expression_bin)+abs(detection_bin-ref$detection_bin)]
  choices<-available[distance==min(distance)][order(gene)]
  chosen<-choices[sample.int(nrow(choices),1L)];used<-c(used,chosen$gene)
  nulls[[length(nulls)+1L]]<-data.table(module=sprintf('null_%03d',m),matched_UPR_gene=g,gene=chosen$gene,expression_bin=chosen$expression_bin,detection_bin=chosen$detection_bin,bin_distance=chosen$distance)
 }
 stopifnot(length(used)==90L,!anyDuplicated(used),!any(used%in%upr))
}
nulls<-rbindlist(nulls);fwrite(nulls,file.path(out,'null100_module_membership.tsv'),sep='\t')
spec<-list(UPR=list(positive=upr,negative=character()))
for(m in unique(nulls$module))spec[[m]]<-list(positive=nulls[module==m,gene],negative=character())
membership<-fread(file.path(cross,'final_axis_membership.tsv'))
for(a in c('apoptosis','proliferation','stress'))spec[[paste0('core_',a)]]<-list(positive=intersect(membership[program==a&direction=='positive_genes',gene],universe),negative=intersect(membership[program==a&direction=='negative_genes',gene],universe))
spec$E2F<-list(positive=intersect(membership[program=='hallmark_E2F_TARGETS',gene],universe),negative=character())
saveRDS(spec,file.path(out,'phase2_axis_spec.rds'))
saveRDS(list(source_path=source_path,universe=universe,original_universe=original_universe,original_upr=original_upr,spec=spec,allgenes=allgenes,raw_gene=raw_gene,source_cells=meta,analysis_genes=eligible$gene,target_genes=target_genes,fixed_external_genes=fixed),file.path(out,'source_preparation.rds'),compress='gzip')
writeLines(c('Metadata-only overlap:454 genes under weak source>=100cells/3GEM-label rule.',
'Prespecified source support:>=3GEMs each>=5targetcells excludes16 of those genes, leaving438.',
'No RNA effect, sign, significance or target-expression gate was used.',
'390 is a separate stricter RPE1>=3GEMs each>=5targetcells sensitivity, not the primary panel.',
'The author RPE1 roster has no QC/doublet flag; finite positive recorded UMI metadata is retained without a fitted classifier.',
'Source cell ordering and gene/GEM annotations were verified against the raw source object.'),file.path(out,'metadata_accounting_notes.txt'))
cat('Source-control preparation and100fixed null memberships complete. RPE1 X remains unread.\n')
