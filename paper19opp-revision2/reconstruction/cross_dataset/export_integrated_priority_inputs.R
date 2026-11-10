#!/usr/bin/env Rscript
# Prespecified secondary integration input, before external response inspection.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
d<-as.data.table(readRDS(file.path(out,'common_cell_scores.rds')))
axes<-c('apoptosis','proliferation','stress')
sd_control<-vapply(axes,function(a)sd(d[dataset=='Replogle'&is_control==TRUE][[paste0('score_',a)]]),numeric(1))
stopifnot(all(is.finite(sd_control)&sd_control>0))
p<-fread(file.path(out,'response_predictions_by_axis_primary.tsv.gz'))[dataset=='Replogle'&model%in%c('pilot_plus_state','pilot_response')&axis%in%axes]
u<-fread(file.path(out,'analysis_units_partition1.tsv.gz'))[dataset=='Replogle',.(unit_id,weight=n_treated_A+n_treated_B,n_treated_A,n_treated_B)]
p<-merge(p,u,by='unit_id');p[,control_sd:=sd_control[axis]]
p[,predicted_z:=prediction/control_sd]
summary<-p[,.(pred_raw=sum(prediction*weight)/sum(weight),pred_z=sum(predicted_z*weight)/sum(weight),n_units=uniqueN(unit_id),covered_treated_cells=sum(weight),pilot_cells=sum(n_treated_A)),by=.(perturbation,model,axis)]
w<-dcast(summary,perturbation+model~axis,value.var=c('pred_raw','pred_z'))
meta<-unique(summary[,.(perturbation,model,n_units,covered_treated_cells,pilot_cells)])
w<-merge(w,meta,by=c('perturbation','model'));setnames(w,'perturbation','gene')
fwrite(w,file.path(out,'integrated_priority_gene_predictions.tsv'),sep='\t')
fwrite(data.table(axis=axes,global_Replogle_control_SD=sd_control),file.path(out,'integrated_priority_control_scales.tsv'),sep='\t')
writeLines(c('Fixed secondary model: source-heldout Replogle pilot_plus_state ridge trained ONLY on Norman and MIX-Seq.',
'All feature scaling and lambda selection use those two training studies, including training-study-only inner tuning.',
'Only primary partition1 is used; no external response outcome is read by this script.',
'Pool predicted matched-control deltas across eligible GEM/target units using treated-cell counts; standardize by harmonized global Replogle control SD.',
'Candidate target program, sign and technical agreement remain fixed by Replogle discovery.',
'Priority is max(fixed discovery sign * predicted target z,0) * max(fixed technical agreement,0)/(1+Euclidean norm of the other two predicted z values).',
'Coverage restriction: at least one metadata-eligible GEM/target unit; do not infer values for uncovered genes.',
'Primary magnitude/composite comparators must be restricted to the identical covered external genes.'),file.path(out,'integrated_priority_lock.txt'))
cat('Integrated gene coverage:',uniqueN(w$gene),'\n')
