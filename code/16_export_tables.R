#!/usr/bin/env Rscript
# Compact, machine-readable core results for direct inspection after reproduction.
suppressPackageStartupMessages(library(data.table))
dir.create("tables",showWarnings=FALSE)
exports <- c(
 task_catalogue="results/task_registry/state_transition_tasks.tsv",
 gene_evidence="results/evidence_validation/gene_level_evidence_status.tsv",
 interactions="results/norman_control/combination_interactions.tsv",
 mixseq_predictions="results/mixseq_control_models/predictions.tsv",
 mixseq_validation_splits="results/mixseq_control_models/split_assignments.tsv",
 combination_predictions="results/grouped_validation/norman_predictions.tsv",
 combination_baseline_predictions="results/extended_baselines/norman_baseline_predictions.tsv",
 multiome_support="results/multiperturb_multiome/multiome_target_audit.tsv",
 norman_growth="results/prediction_validation/functional_condition_mapping.tsv",
 lineage_associations="results/lineage_primary/state_associations.tsv",
 lineage_transfer="results/lineage_primary/leave_dataset_transfer.tsv"
)
for (nm in names(exports)) {
  x <- fread(exports[[nm]])
  if(nm=="gene_evidence") {
    q<-fread("results/multiome_validation/gene_criterion_four_state.tsv")
    stopifnot(all(merge(x[,.(gene,evidence_class)],q[,.(gene,evidence_class)],by="gene")[,evidence_class.x==evidence_class.y]))
    q<-q[,c("gene",setdiff(names(q),names(x))),with=FALSE]
    x<-merge(x,q,by="gene",sort=FALSE)
    old_status<-intersect(c("causal_movement_status","two_context_status","combination_status","specificity_status"),names(x));x[,(old_status):=NULL]
    changes<-fread("results/multiome_validation/nt_required_evidence_changed_genes.tsv")
    x[,evidence_class_with_NT_requirement:=evidence_class]
    x[changes,on="gene",evidence_class_with_NT_requirement:=i.sensitivity_class]
    x[,contexts_with_NT_requirement:=independent_contexts_with_causal_support-as.integer(multiperturb_eligible)]
    stopifnot(sum(x$evidence_class_with_NT_requirement=="causal_regulator")==389L)
  }
  if(nm=="task_catalogue") {
    nt<-fread("results/multiome_validation/non_targeting_target_comparison.tsv")[nt_construct_fdr<.05,gene]
    x[,eligible_if_multiome_NT_required:=status=="eligible"&(dataset_id!="yan_2024_multiperturbseq"|perturbation_identity%in%nt)]
    stopifnot(sum(x$eligible_if_multiome_NT_required)==705L)
  }
  dest <- gzfile(file.path("tables",paste0(nm,".tsv.gz")),"wt",compression=9)
  write.table(x,dest,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
  close(dest)
}
stopifnot(sum(fread(exports[["task_catalogue"]])$status=="eligible")==715,
          nrow(fread(exports[["gene_evidence"]]))==9894,
          nrow(fread(exports[["interactions"]]))==131)
