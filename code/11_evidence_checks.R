#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(TRUE); base<-"."; out<-"results/evidence_validation"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
e<-fread(file.path(base,"results/control_atlas/control_node_evidence.tsv"))
t<-fread(file.path(base,"results/task_registry/state_transition_tasks.tsv"))
n<-fread(file.path(base,"results/norman_control/condition_index.tsv"))
m<-fread(file.path(base,"results/multiperturb_multiome/target_reproducibility.tsv"))
replogle_genes<-unique(t[dataset_id=="replogle_2022_k562_essential",perturbation_identity])
norman_genes<-unique(n[condition_type=="single" & !is.na(single_gene),single_gene])
multi_genes<-unique(m$guide_gene)
e[,`:=`(replogle_assayed=gene%in%replogle_genes,norman_assayed=gene%in%norman_genes,multiome_assayed=gene%in%multi_genes)]
e[,supporting_datasets:=vapply(seq_len(.N),function(i)paste(c("Replogle_essential_CRISPRi_K562","Norman_CRISPRa_K562","MultiPerturb_CRISPRi_ATRT")[c(replogle_eligible[i],norman_single_eligible[i],multiperturb_eligible[i])],collapse=";"),character(1))]
e[,assayed_contexts:=vapply(seq_len(.N),function(i) paste(c("K562","ATRT_human")[c(replogle_assayed[i]|norman_assayed[i],multiome_assayed[i])],collapse=";"),character(1))]
e[,causal_contexts:=vapply(seq_len(.N),function(i) paste(c("K562","ATRT_human")[c(replogle_eligible[i]|norman_single_eligible[i],multiperturb_eligible[i])],collapse=";"),character(1))]
e[,perturbation_modality_note:="CRISPRi and CRISPRa are counted as separate source evidence; their effect signs are not assumed equivalent"]
e[,causal_movement_status:=fifelse(causal_regulator,"met",fifelse(replogle_assayed|norman_assayed|multiome_assayed,"not_met_under_QC_rules","not_evaluable_no_cell_level_target_task"))]
e[,two_context_status:=fifelse(independent_contexts_with_causal_support>=2,"met",fifelse((replogle_assayed|norman_assayed)&multiome_assayed,"not_met_under_QC_rules","not_evaluable_no_second_context"))]
e[,combination_status:=fifelse(observed_combinations==0,"not_evaluable_no_pair_measurement",fifelse(combination_consistent,"met","not_met"))]
e[,specificity_status:=fifelse(!(replogle_eligible|norman_single_eligible),"not_evaluable_no_eligible_K562_response",fifelse(path_specific,"met","not_met"))]
e[,all_candidate_evidence_available:=two_context_status!="not_evaluable_no_second_context" & observed_combinations>0 & (replogle_eligible|norman_single_eligible)]
e[,assayed_datasets:=vapply(seq_len(.N),function(i)paste(c("Replogle_essential","Norman","MultiPerturb","Replogle_genomewide_breadth")[c(replogle_assayed[i],norman_assayed[i],multiome_assayed[i],genomewide_screened[i])],collapse=";"),character(1))]
e[,next_evidence_needed:=vapply(seq_len(.N),function(i) paste(c(if(!causal_regulator[i])"quality-supported target perturbation with independent construct or experiment",if(two_context_status[i]!="met")"eligible target response in a second independent biological context",if(combination_status[i]!="met")"measured combination response consistent with a prespecified target",if(specificity_status[i]!="met")"target-specific response with quantified non-target deviation", "held-out fixed-phenotype validation and durability/washout test before validated-control claim"),collapse=";"),character(1))]
setorder(e,gene)
fwrite(e,file.path(out,"gene_level_evidence_status.tsv"),sep="\t",na="NA")
fwrite(e[genomewide_screened==FALSE,.(gene,assayed_contexts,supporting_datasets)],file.path(out,"additional_genes_outside_breadth_screen.tsv"),sep="\t")
status_counts<-rbindlist(lapply(c("causal_movement_status","two_context_status","combination_status","specificity_status"),function(z) e[,.(criterion=z,status=get(z))][,.N,by=.(criterion,status)]))
fwrite(status_counts,file.path(out,"evidence_criterion_counts.tsv"),sep="\t")
sens<-rbindlist(lapply(c(1,2),function(ctx)rbindlist(lapply(c(.25,.5,.75),function(q)rbindlist(lapply(c(.3,.5,.7),function(cos){
 cut<-quantile(e$off_target_deviation[e$replogle_eligible],q,na.rm=TRUE)
 spec<-(e$replogle_eligible & e$off_target_deviation<=cut)|(e$norman_single_eligible&e$norman_single_reproducibility>=cos)
 for_comb<-list(primary=e$combination_consistent,at_least_one_positive=e$synergistic_combinations>=1,nonantagonistic_positive=e$synergistic_combinations>=1&e$antagonistic_combinations==0)
 rbindlist(lapply(names(for_comb),function(rule) data.table(minimum_contexts=ctx,replogle_offtarget_quantile=q,norman_reproducibility_cosine=cos,combination_rule=rule,n_candidates=sum(e$causal_regulator&e$independent_contexts_with_causal_support>=ctx&for_comb[[rule]]&spec,na.rm=TRUE))))
}))))))
fwrite(sens,file.path(out,"candidate_criteria_sensitivity.tsv"),sep="\t")
summary<-data.table(metric=c("union_genes","breadth_screen_genes","extra_genes","causal_regulators","state_markers","insufficient_causal","second_context_assayed","second_context_causal_supported","combination_observed","all_candidate_evidence_available","candidates_primary"),value=c(nrow(e),sum(e$genomewide_screened),sum(!e$genomewide_screened),sum(e$causal_regulator),sum(e$evidence_class=="state_marker"),sum(e$evidence_class=="insufficient_causal_evidence"),sum((e$replogle_assayed|e$norman_assayed)&e$multiome_assayed),sum(e$independent_contexts_with_causal_support>=2),sum(e$observed_combinations>0),sum(e$all_candidate_evidence_available),sum(e$validated_control_node_candidate)))
fwrite(summary,file.path(out,"evidence_universe_summary.tsv"),sep="\t")
print(summary);print(status_counts);print(sens[,.(minimum=min(n_candidates),maximum=max(n_candidates)),by=minimum_contexts]);print(e[independent_contexts_with_causal_support>=2,.(gene,combination_status,specificity_status)])
