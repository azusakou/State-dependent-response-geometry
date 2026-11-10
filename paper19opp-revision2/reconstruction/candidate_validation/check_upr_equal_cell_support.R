#!/usr/bin/env Rscript
# Sampling and guide-assignment sensitivities of the exploratory UPR association.
suppressPackageStartupMessages(library(data.table));a<-commandArgs(TRUE);root<-if(length(a))normalizePath(a[1])else normalizePath('.');out<-file.path(root,'paper19opp-revision2/reconstruction/candidate_validation');set.seed(202611063L);B<-500L
rep<-as.data.table(readRDS(file.path(root,'paper19opp-revision2/reconstruction/cross_dataset/harmonized_score_part_1.rds')));rep[,gene:=as.character(gene)];u<-'score_hallmark_UNFOLDED_PROTEIN_RESPONSE';axiscols<-paste0('score_',c('apoptosis','proliferation','stress','hallmark_UNFOLDED_PROTEIN_RESPONSE','hallmark_E2F_TARGETS','hallmark_G2M_CHECKPOINT','hallmark_P53_PATHWAY','hallmark_TNFA_SIGNALING_VIA_NFKB','hallmark_INTERFERON_GAMMA_RESPONSE','hallmark_MTORC1_SIGNALING','hallmark_OXIDATIVE_PHOSPHORYLATION','hallmark_HYPOXIA'))
rr<-list();summary<-list();counts<-list()
run<-function(z,st,cohort){for(ctx in unique(z$context))for(mode in c('original_assignments','single_guide_assignments')){
 if(mode=='single_guide_assignments'&&!('simple_single_guide'%in%names(z)))next
 d<-z[context==ctx];if(mode=='single_guide_assignments')d<-d[simple_single_guide==TRUE]
 for(m in c(30L,15L)){
  cnt<-d[is_control==FALSE,.N,by=gene];gg<-intersect(cohort,cnt[N>=m,gene]);gg<-intersect(gg,rep[gene!='non-targeting',.N,by=gene][N>=m,gene]);n<-length(gg);counts[[length(counts)+1L]]<<-data.table(study=st,context=ctx,assignment=mode,cells_per_target=m,genes=n,control_cells=d[is_control==TRUE,.N]);if(n<10||d[is_control==TRUE,.N]<m)next
  sr<-rep[gene%in%gg];vt<-d[gene%in%gg&is_control==FALSE];sindex<-split(seq_len(nrow(sr)),sr$gene);vindex<-split(seq_len(nrow(vt)),vt$gene)
  for(b in seq_len(B)){
   sm<-vm<-matrix(NA_real_,n,12);for(i in seq_along(gg)){is<-sindex[[gg[i]]];iv<-vindex[[gg[i]]];is<-is[sample.int(length(is),m,replace=FALSE)];iv<-iv[sample.int(length(iv),m,replace=FALSE)];sm[i,]<-colMeans(as.matrix(sr[is,..axiscols]));vm[i,]<-colMeans(as.matrix(vt[iv,..axiscols]))}
   rh<-vapply(seq_len(12),function(j)cor(sm[,j],vm[,j],method='spearman'),numeric(1));rr[[length(rr)+1L]]<<-data.table(study=st,context=ctx,assignment=mode,cells_per_target=m,genes=n,draw=b,axis=sub('^score_','',axiscols),spearman=rh)
  }
 }
}}
adam<-fread(file.path(out,'adamson_cell_scores_harmonized.tsv.gz'));fr<-fread(file.path(out,'frangieh_cell_scores_harmonized.tsv.gz'));run(adam,'Adamson2016',fread(file.path(out,'harmonized_adamson_candidates.tsv'))$gene);run(fr,'Frangieh2021',fread(file.path(out,'harmonized_frangieh_candidates.tsv'))$gene)
r<-rbindlist(rr);fwrite(r,file.path(out,'exploratory_equal_cell_axis_correlations.tsv.gz'),sep='\t');s<-r[,.(draws=.N,median_spearman=median(spearman),sampling_percentile_low=quantile(spearman,.025),sampling_percentile_high=quantile(spearman,.975),positive_fraction=mean(spearman>0)),by=.(study,context,assignment,cells_per_target,genes,axis)];fwrite(s,file.path(out,'exploratory_equal_cell_axis_summary.tsv'),sep='\t');fwrite(rbindlist(counts),file.path(out,'exploratory_equal_cell_metadata_counts.tsv'),sep='\t');print(s[axis=='hallmark_UNFOLDED_PROTEIN_RESPONSE'])
