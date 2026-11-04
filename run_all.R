#!/usr/bin/env Rscript
# Numerical reproduction from compact, immutable source-derived inputs.
a<-commandArgs(TRUE);quick<-'--quick'%in%a
needed<-c('data.table','cluster','glmnet','nnet','igraph','xgboost','BiocNeighbors','R.utils','Matrix')
missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Install required packages: ',paste(missing,collapse=', '))
if(!file.exists('inputs/auxiliary.rds'))stop('Run this command from the package root.')
dir.create('results',showWarnings=FALSE)
for(stage in c('01_mixseq_probabilities.R','02_norman_pairs.R','03_multiome.R','04_task_registry.R','05_mixseq_models.R','06_grouped_models.R','07_extended_baselines.R','08_evidence_graphs.R','09_prediction_checks.R','11_evidence_checks.R','12_norman_checks.R','13_lineage_additional_checks.R','14_prediction_additional_checks.R','15_multiome_additional_checks.R','16_export_tables.R')){
 cat('\nRunning',stage,'\n');source(file.path('code',stage),local=new.env(parent=globalenv()))
}
if(!quick){cat('\nRunning full clone-neighbourhood sensitivities and 200-draw bootstrap\n');source('code/10_lineage.R',local=new.env(parent=globalenv()))}
# Append new learned-representation checks after every original stage, including
# the full lineage bootstrap, so their RNG cannot change any published replay.
cat('\nRunning fully nested scVI prediction with two prespecified NN budgets\n')
nested<-new.env(parent=globalenv());source('code/17_nested_scvi_models.R',local=nested)
n600<-nested$run_nested_scvi_models(getwd(),'inputs/mixseq_nested_scvi.rds','results/nested_scvi_600',nn_maxit=600L)
n6000<-nested$run_nested_scvi_models(getwd(),'inputs/mixseq_nested_scvi.rds','results/nested_scvi_6000',nn_maxit=6000L)
nested$compare_nested_scvi_nn_budgets('results/nested_scvi_600','results/nested_scvi_6000','results/nested_scvi_nn6000',n600,n6000)
nested$export_nested_scvi_tables(getwd(),'results/nested_scvi_600','results/nested_scvi_6000','results/nested_scvi_nn6000','results/nested_scvi')
writeLines(capture.output(sessionInfo()),'results/sessionInfo.txt')
source('code/verify.R',local=new.env(parent=globalenv()))
cat('\nCompleted requested reproduction stages.\n')
