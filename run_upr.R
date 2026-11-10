#!/usr/bin/env Rscript
# Reproduce frozen UPR metrics from compact, source-derived inputs.
suppressPackageStartupMessages(library(data.table))
stopifnot(file.exists('inputs/upr/HepG2_gene_effects.tsv'))
if(!requireNamespace('digest',quietly=TRUE))stop('Install digest.')
source('code/18_upr_holdout.R',local=new.env(parent=globalenv()))
source('code/19_upr_gem_sensitivity.R',local=new.env(parent=globalenv()))
for(pair in list(c('primary_Holm5.tsv','upr_primary_Holm5.tsv'),c('all_ordering_results.tsv','upr_ordering.tsv'),c('all_fixed_baseline_metrics.tsv','upr_baseline_metrics.tsv'),c('shared_GEM_intervals.tsv','upr_GEM_intervals.tsv'))){
 actual<-fread(file.path('results/upr',pair[1]));expected<-fread(file.path('reference',pair[2]))
 stopifnot(identical(names(actual),names(expected)),nrow(actual)==nrow(expected))
 for(nm in names(actual)){
  x<-actual[[nm]];y<-expected[[nm]]
  if(is.numeric(x)&&is.numeric(y)){
   stopifnot(identical(is.na(x),is.na(y)),all(abs(x[is.finite(x)]-y[is.finite(y)])<=1e-10))
  }else stopifnot(identical(as.character(x),as.character(y)))
 }
}
cat('All frozen UPR numerical references reproduced.\n')
