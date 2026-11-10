#!/usr/bin/env Rscript
# Metadata-only reconstruction; no X values are accessed in this stage.
project<-normalizePath(commandArgs(TRUE)[1]);out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'))
old<-file.path(project,'paper19opp-revision2/reconstruction/cross_dataset/phase2_upr_rpe1');p<-readRDS(file.path(old,'source_preparation.rds'))
paths<-c(K562=p$source_path,Adamson='/Volumes/wang/paper4rfv/data/processed/h5ad/core_unpacked/adamson_2016_upr__GSM2406681_10X010.h5ad',Frangieh='/Volumes/wang/paper13/raw/perturbation/SCP1064/FrangiehIzar2021_RNA.h5ad',RPE1='/Volumes/wang/paper6yhn/raw/replogle_rpe1/rpe1_raw_singlecell_01.h5ad')
maps<-sets<-list()
for(nm in names(paths)){
 f<-hdf5r::H5File$new(paths[nm],'r');fld<-if(nm=='Frangieh')'gene_symbol'else if(nm%in%c('K562','RPE1'))'gene_name'else'_index';if(!fld%in%names(f[['var']]))fld<-'gene'
 g<-as.character(annotation(f,'var',fld));sets[[nm]]<-unique(g);maps[[nm]]<-data.table(study=nm,raw_feature_row=seq_along(g),gene=g);f$close_all()
}
fwrite(rbindlist(maps),file.path(out,'four_study_feature_mapping.tsv'),sep='\t');common4<-sort(Reduce(intersect,sets));fwrite(data.table(gene=common4),file.path(out,'four_study_metadata_intersection.tsv'),sep='\t')
primary<-fread(file.path(old,'phase2_gene_universe.tsv'))$gene;stopifnot(length(primary)==6780L,all(primary%in%common4));fwrite(data.table(gene=primary),file.path(out,'primary_gene_universe.tsv'),sep='\t')
saveRDS(list(paths=paths,rawmaps=maps,universe=primary,source_cells=p$source_cells),file.path(out,'metadata_preparation.rds'))
