#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix);library(digest)})
project<-normalizePath(commandArgs(TRUE)[1]);out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');old<-file.path(project,'paper19opp-revision2/reconstruction/cross_dataset/phase2_upr_rpe1')
locks<-fread(file.path(out,'membership_lock_manifest.tsv'));stopifnot(all(vapply(file.path(out,locks$file),digest,character(1),file=TRUE,algo='sha256')==locks$sha256))
membership<-fread(file.path(out,'locked_memberships.tsv'));membership[,measured:=as.logical(measured)]
defs<-fread(file.path(out,'locked_axis_definitions.tsv'));universe<-fread(file.path(out,'primary_gene_universe.tsv'))$gene
spec<-setNames(lapply(defs$axis,function(a)list(positive=membership[axis==a&measured==TRUE,assay_gene],negative=character())),defs$axis)
oldspec<-readRDS(file.path(old,'phase2_axis_spec.rds'));spec$Hallmark_UPR90<-oldspec$UPR
for(a in c('core_apoptosis','core_proliferation','core_stress','E2F'))spec[[a]]<-oldspec[[a]]
stopifnot(length(spec)==21L,length(spec$Hallmark_UPR90$positive)==90L,all(unlist(lapply(spec,function(a)c(a$positive,a$negative)))%in%universe))
pool<-fread(file.path(old,'source_NTC_matching_bins.tsv'));stopifnot(nrow(pool)==6780L,setequal(pool$gene,universe),all(pool$source_NTC_cells==9775L))
fwrite(pool,file.path(out,'source_NTC_expression_detection_bins.tsv'),sep='\t')
amap<-fread(file.path(out,'assay_HGNC_mapping.tsv'));hfull<-fread(file.path(out,'original_Hallmark_UPR_full.tsv'));hallmark_exclude<-amap[canonical_gene%in%hfull$canonical_gene&canonical_gene!='',assay_gene]
stopifnot(all(spec$Hallmark_UPR90$positive%in%hallmark_exclude))
families<-c('ER_folding_full','ERAD_full','ER_targeting_full','ISR_full','Reactome_UPR_full','GO_ER_UPR_full','Reactome_UPR_no_Hallmark','GO_ER_UPR_no_Hallmark')
draws<-list()
for(fi in seq_along(families)){
 family<-families[fi];members<-sort(spec[[family]]$positive)
 for(m in 1:100){
  set.seed(20261106L+fi*1000L+m);used<-character()
  for(g in members){
   ref<-pool[gene==g];available<-pool[!gene%in%c(members,hallmark_exclude,used)]
   available[,distance:=abs(expression_bin-ref$expression_bin)+abs(detection_bin-ref$detection_bin)]
   choices<-available[distance==min(distance)][order(gene)];stopifnot(nrow(choices)>0L)
   chosen<-choices[sample.int(nrow(choices),1L)];used<-c(used,chosen$gene)
   draws[[length(draws)+1L]]<-data.table(family=family,module=sprintf('null__%s__%03d',family,m),seed=20261106L+fi*1000L+m,matched_gene=g,gene=chosen$gene,expression_bin=chosen$expression_bin,detection_bin=chosen$detection_bin,bin_distance=chosen$distance)
  }
  stopifnot(length(used)==length(members),!anyDuplicated(used),!any(used%in%c(members,hallmark_exclude)))
 }
}
nulls<-rbindlist(draws);fwrite(nulls,file.path(out,'null100_per_readout_memberships.tsv.gz'),sep='\t')
for(m in unique(nulls$module))spec[[m]]<-list(positive=nulls[module==m,gene],negative=character())
summary<-nulls[,.(random_sets=uniqueN(module),readout_members=uniqueN(matched_gene),member_draws=.N,unique_sampled_genes=uniqueN(gene),nonexact_bin_draws=sum(bin_distance>0),max_bin_distance=max(bin_distance)),by=family]
fwrite(summary,file.path(out,'null_matching_accounting.tsv'),sep='\t')
fwrite(data.table(assay_gene=hallmark_exclude),file.path(out,'Hallmark_disjoint_measured_exclusion.tsv'),sep='\t')
saveRDS(spec,file.path(out,'scoring_spec.rds'))
fixed438<-fread(file.path(old,'phase2_eligible_genes.tsv'))$gene
adam<-fread(file.path(project,'paper19opp-revision2/reconstruction/candidate_validation/harmonized_adamson_candidates.tsv'))$gene
frang<-fread(file.path(project,'paper19opp-revision2/reconstruction/candidate_validation/harmonized_frangieh_candidates.tsv'))$gene
roster<-rbindlist(list(data.table(study='K562',gene=sort(unique(c(fixed438,adam,frang)))),data.table(study='RPE1',gene=fixed438),data.table(study='Adamson',gene=adam),data.table(study='Frangieh',gene=frang)))
fwrite(roster,file.path(out,'fixed_study_gene_rosters.tsv'),sep='\t')
writeLines(c(paste('Frozen UTC',format(Sys.time(),tz='UTC',usetz=TRUE)),
'Source-only controls are the same9775K562non-targeting cells and fixed6780CP10k universe used in the preserved source-control statistics. Their exact source-only expression/detection summaries and5x5rank-based bins are reused without reference to validation expression.',
'Eight families each have100 matched random modules. Module size equals its measured real-readout size; samples are without replacement within a module. The tested full/readout members and all full-Hallmark canonical-overlap assay genes are excluded from its pool. Each member matches source-control expression/detection bins; nearest Manhattan bin distance is used only if the exact bin has no available gene. All distances and seeds are retained. No outcomebased replacement or selection occurs.',
'Real signatures are scored in exactL-minus-target primary and fixedL6780 sensitivity modes. Nulls are fixedL6780 only, with the identical target-numerator correction as the paired real readout. Exact-denominator-primary null specificity is not estimated by these nulls.'),file.path(out,'null_matching_protocol.txt'))
files<-c('membership_lock_manifest.tsv','null100_per_readout_memberships.tsv.gz','null_matching_accounting.tsv','source_NTC_expression_detection_bins.tsv','scoring_spec.rds','fixed_study_gene_rosters.tsv','null_matching_protocol.txt','prepare_nulls.R')
manifest<-data.table(file=files,sha256=vapply(file.path(out,files),digest,character(1),file=TRUE,algo='sha256'))
fwrite(manifest,file.path(out,'scoring_lock_manifest.tsv'),sep='\t')
cat('Scoring/null lock SHA256',digest(file.path(out,'scoring_lock_manifest.tsv'),file=TRUE,algo='sha256'),'\n');print(summary);print(roster[,.N,by=study]);writeLines(capture.output(sessionInfo()),file.path(out,'preparation_sessionInfo.txt'))
