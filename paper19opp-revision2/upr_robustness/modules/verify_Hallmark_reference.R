#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(args[1]);study_name<-args[2]
suppressPackageStartupMessages({library(data.table);library(digest)})
out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');old<-file.path(project,'paper19opp-revision2/reconstruction/cross_dataset/phase2_upr_rpe1');dest<-file.path(out,'scored')
new<-fread(file.path(dest,paste0(study_name,'_gene_effects.tsv.gz')))[axis=='Hallmark_UPR90']
cache<-readRDS(file.path(dest,paste0(study_name,'_cell_cache.rds')))
if(study_name=='K562'){
 fixed<-fread(file.path(old,'source_rank_inputs.tsv'))[,.(gene,reference_delta_z=source_UPR_deltaZ,reference_n=source_cells)]
 exact<-fread(file.path(old,'source_UPR_exact_denominator_gene_effects.tsv'))[,.(gene,reference_delta_z=deltaZ,reference_n=source_cells)]
 raw<-readRDS(file.path(old,'original91_source_raw_cache.rds'));ids<-match(cache$metadata$cell_id,raw$cells$cell_id);stopifnot(!anyNA(ids));Lold<-raw$cells$phase2_library[ids]
}else if(study_name=='RPE1'){
 fixed<-fread(file.path(project,'paper19opp-revision2/reconstruction/upr_rpe1_validation/RPE1_gene_effects_all105.tsv.gz'))[axis=='UPR',.(gene,reference_delta_z=deltaZ,reference_n=source_cells)]
 exact<-fread(file.path(old,'rpe1_UPR_exact_denominator_gene_effects.tsv'))[,.(gene,reference_delta_z=deltaZ,reference_n=source_cells)]
 raw<-readRDS(file.path(old,'RPE1_phase2_cache.rds'));ids<-match(cache$metadata$cell_id,raw$cell_metadata$cell_id);stopifnot(!anyNA(ids));Lold<-raw$cell_metadata$common_library[ids]
}else stop('Reference comparison is defined for original438 source/RPE panels only')
qa<-list()
for(mode_name in c('fixed_L6780','exact_Lminus_target')){
 ref<-if(mode_name=='fixed_L6780')fixed else exact
 z<-merge(new[mode==mode_name],ref,by='gene');stopifnot(nrow(z)==438L,all(z$n_target==z$reference_n),all(is.finite(z$delta_z)))
 qa[[mode_name]]<-data.table(study=study_name,mode=mode_name,reference_genes=nrow(z),max_abs_effect_difference=max(abs(z$delta_z-z$reference_delta_z)),max_abs_common_library_difference=max(abs(cache$metadata$common_library-Lold)),same_target_counts=all(z$n_target==z$reference_n),passed=max(abs(z$delta_z-z$reference_delta_z))<1e-8&&max(abs(cache$metadata$common_library-Lold))==0)
}
qa<-rbindlist(qa);fwrite(qa,file.path(dest,paste0(study_name,'_Hallmark_reference_QA.tsv')),sep='\t');stopifnot(all(qa$passed));print(qa)
