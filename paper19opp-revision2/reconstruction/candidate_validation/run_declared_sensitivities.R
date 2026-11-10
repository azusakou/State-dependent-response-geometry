#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(TRUE);out<-if(length(args))args[1]else'paper19opp-revision2/reconstruction/candidate_validation'
s<-readLines(file.path(out,'analyze_locked_validation.R'));cut<-grep('^# Required primary external tables',s)[1];eval(parse(text=paste(s[seq_len(cut-1)],collapse='\n')))
adam<-fread(file.path(out,'adamson_cell_scores_harmonized.tsv.gz'));ae<-fread(file.path(out,'adamson_control_target_expression_harmonized.tsv'));fr<-fread(file.path(out,'frangieh_cell_scores_harmonized.tsv.gz'));fe<-fread(file.path(out,'frangieh_control_target_expression_harmonized.tsv'));frgenes<-fread(file.path(out,'harmonized_frangieh_candidates.tsv'))$gene
# Planned comparator availability is metadata-defined; no rank or endpoint is refitted.
fitnessgenes<-rank[is.finite(depmap_generic_essentiality),gene]
process(adam,ae,'Adamson2016 essentiality-matched','is_control',100L,fitnessgenes)
process(fr,fe,'Frangieh2021 essentiality-matched','is_control',30L,intersect(frgenes,fitnessgenes))
expressiongenes<-rank[is.finite(discovery_control_expression),gene]
process(fr,fe,'Frangieh2021 expression-matched','is_control',30L,intersect(frgenes,expressiongenes))
# Reference and normalization sensitivities were specified before response evaluation.
process(adam,ae,'Adamson2016 alternative-control sensitivity','is_alternative_control',100L)
fa<-copy(adam);for(ax in axes)fa[,(paste0('score_',ax)):=get(paste0('full_library_score_',ax))]
process(fa,ae,'Adamson2016 full-library sensitivity','is_control',100L)
ff<-copy(fr);for(ax in axes)ff[,(paste0('score_',ax)):=get(paste0('score_full_library_',ax))]
process(ff,fe,'Frangieh2021 full-library sensitivity','is_control',30L,frgenes)
for(gate in c(20L,50L,100L))process(fr,fe,paste0('Frangieh2021 count-gate-',gate,' sensitivity'),'is_control',gate)
for(z in c('metric_rows','gene_rows','choice_rows','equal_rows','contrast_rows','null_rows'))fwrite(rbindlist(get(z),fill=TRUE),file.path(out,paste0(sub('_rows$','',z),'_sensitivity.tsv.gz')),sep='\t')
fwrite(data.table(issue=c('DepMap primary Adamson coverage','DepMap primary Frangieh coverage','Discovery baseline-expression Frangieh coverage'),covered=c(sum(fread(file.path(out,'harmonized_adamson_candidates.tsv'))$gene%in%fitnessgenes),sum(frgenes%in%fitnessgenes),sum(frgenes%in%expressiongenes)),total=c(30,19,19)),file.path(out,'comparator_availability.tsv'),sep='\t')
