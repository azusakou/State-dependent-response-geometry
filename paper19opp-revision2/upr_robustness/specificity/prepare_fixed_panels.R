#!/usr/bin/env Rscript
# Metadata-only projection of the existing target cohorts; no new response values.
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);project<-if(length(a))a[1]else'paper19opp-revision2'
d<-fread(file.path(project,'reconstruction/candidate_validation/exploratory_all_gene_axis_effects.tsv.gz'),select=c('study','context','gene','axis'))[axis=='hallmark_UNFOLDED_PROTEIN_RESPONSE']
d[,comparison:=fifelse(study=='Adamson2016','Adamson',fifelse(context=='Control','Frangieh_Control',fifelse(context=='IFNγ','Frangieh_IFNg','Frangieh_Coculture')))];d[,validation_study:=fifelse(study=='Adamson2016','Adamson','Frangieh')];d[,validation_context:=context]
r<-fread(file.path(project,'reconstruction/cross_dataset/phase2_upr_rpe1/source_rank_inputs.tsv'),select='gene');r[,`:=`(comparison='RPE1',validation_study='RPE1',validation_context='RPE1')]
z<-rbindlist(list(d[,.(comparison,validation_study,validation_context,gene)],r),use.names=TRUE);setorder(z,comparison,gene)
stopifnot(nrow(z)==523L,!anyDuplicated(z[,.(comparison,gene)]),identical(sort(z[,.N,by=comparison]$N),c(17L,19L,19L,30L,438L)))
fwrite(z,file.path(project,'upr_robustness/specificity/fixed_target_panels.tsv'),sep='\t');print(z[,.N,by=comparison])
