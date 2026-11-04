#!/usr/bin/env Rscript

# Atlas contract
# Core conclusion: target basins differ in control degeneracy, while only a subset
# of causal perturbations satisfy the full control-node evidence ladder.
# Sequential policies remain inside observed time/intervention support.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})
source("code/common.R")

set.seed(1906)
out_dir <- "results/control_atlas"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

tasks <- fread("results/task_registry/state_transition_tasks.tsv")
norman_interaction <- fread(
  "results/norman_control/combination_interactions.tsv"
)
norman_response <- as.data.table(readRDS("results/norman_control/condition_response_delta.rds"))

# Control degeneracy and substitutability among combinations reaching the same
# independently learned response basin.
norman_tasks <- tasks[
  dataset_id == "norman_2019_combinatorial" &
    intervention_size == 2L & status == "eligible",
  .(combination = perturbation_identity, target_state_id, control_energy,
    target_state_distance)
]
norman_tasks <- merge(norman_tasks, norman_interaction,
  by = "combination", all.x = TRUE
)
norman_tasks[, control_energy := control_energy.x]
feature_cols <- setdiff(names(norman_response), "condition_id")
response_mat <- as.matrix(norman_response[, ..feature_cols])
rownames(response_mat) <- norman_response$condition_id

calculate_response_path_diversity <- function(cosine_threshold) {
  rbindlist(lapply(
    split(norman_tasks, norman_tasks$target_state_id),
    function(dt) {
    ids <- dt$combination
    mat <- response_mat[ids, , drop = FALSE]
    norms <- sqrt(rowSums(mat^2))
    cosine <- tcrossprod(mat) / pmax(outer(norms, norms), 1e-12)
    diag(cosine) <- 1
    graph <- graph_from_adjacency_matrix(
      cosine >= cosine_threshold, mode = "undirected", diag = FALSE
    )
    components <- components(graph)$membership
    weight <- dt$target_state_distance / pmax(dt$control_energy, 1e-8)
    component_weight <- rowsum(weight, components, reorder = FALSE)[, 1]
    p <- component_weight / sum(component_weight)
    effective_paths <- exp(-sum(p * log(p)))

    genes <- strsplit(ids, "_", fixed = TRUE)
    pair_index <- which(upper.tri(cosine), arr.ind = TRUE)
    disjoint <- vapply(seq_len(nrow(pair_index)), function(i) {
      length(intersect(
        genes[[pair_index[i, 1]]], genes[[pair_index[i, 2]]]
      )) == 0L
    }, logical(1))
    disjoint_cos <- cosine[pair_index][disjoint]
    substitutable <- if (length(disjoint_cos)) {
      mean(disjoint_cos >= cosine_threshold)
    } else {
      NA_real_
    }
    data.table(
      target_state_id = dt$target_state_id[1],
      observed_control_sets = nrow(dt),
      non_equivalent_path_components = uniqueN(components),
      control_degeneracy = effective_paths,
      pathway_substitutability = substitutable,
      median_control_energy = median(dt$control_energy, na.rm = TRUE),
      fragile_control_axis = effective_paths < 2,
      redundant_control_landscape = effective_paths >= 4
    )
    }
  ))[, cosine_threshold := cosine_threshold]
}
degeneracy_sensitivity <- rbindlist(lapply(
  c(0.70, 0.80, 0.90), calculate_response_path_diversity
))
degeneracy <- degeneracy_sensitivity[cosine_threshold == 0.80]
fwrite(degeneracy, file.path(out_dir, "control_degeneracy.tsv"), sep = "\t")
fwrite(
  degeneracy_sensitivity,
  file.path(out_dir, "control_degeneracy_sensitivity.tsv"),
  sep = "\t"
)

# Evidence ladder for genetic control-node candidates.
rep_tasks <- tasks[
  dataset_id == "replogle_2022_k562_essential",
  .(
    gene = perturbation_identity,
    replogle_eligible = status == "eligible",
    target_state_id,
    target_state_distance,
    control_success_rate,
    off_target_deviation,
    control_energy,
    target_reproducibility
  )
]

norman_index <- fread("results/norman_control/condition_index.tsv")
norman_qc <- fread("results/norman_control/condition_qc.tsv")
norman_single <- merge(norman_index, norman_qc,
  by = c("condition_id", "raw_guide_identity", "n_cells", "condition_type")
)[condition_type == "single" & !is.na(single_gene), .(
  norman_single_conditions = .N,
  norman_single_cells = sum(n_cells),
  norman_single_movement = weighted.mean(response_norm, n_cells),
  norman_single_reproducibility = weighted.mean(
    median_replicate_cosine, n_cells
  ),
  norman_single_eligible = any(passes_task_qc & n_cells >= 100L)
), by = .(gene = single_gene)]

combo_gene <- rbindlist(list(
  norman_interaction[, .(
    gene = gene_a,
    interaction_class,
    target_movement,
    off_target_interaction
  )],
  norman_interaction[, .(
    gene = gene_b,
    interaction_class,
    target_movement,
    off_target_interaction
  )]
))[, .(
  observed_combinations = .N,
  synergistic_combinations = sum(interaction_class == "synergistic"),
  antagonistic_combinations = sum(interaction_class == "antagonistic"),
  redundant_combinations = sum(interaction_class == "redundant"),
  median_combination_movement = median(target_movement, na.rm = TRUE),
  median_combination_off_target = median(
    off_target_interaction, na.rm = TRUE
  )
), by = gene]

multi <- fread("results/multiperturb_multiome/target_reproducibility.tsv")
multi_gene <- multi[, .(
  multiperturb_contexts = .N,
  multiperturb_eligible = any(passes_task_qc & n_cells >= 100L),
  multiperturb_rna_movement = safe_max(rna_movement),
  multiperturb_atac_movement = safe_max(atac_movement),
  multiperturb_reproducibility = safe_max(pmin(
    median_rna_guide_cosine, median_atac_guide_cosine
  ))
), by = .(gene = guide_gene)]

# Genome-wide pseudobulk evidence is included as breadth/association evidence,
# not as a state-transition or causal-control task.
gw <- as.data.table(readRDS("inputs/auxiliary.rds")$gw)
gw[, gene := vapply(
  strsplit(gene_transcript, "_", fixed = TRUE),
  function(x) if (length(x) >= 2L) x[[2]] else NA_character_,
  character(1)
)]
gw_gene <- gw[!core_control & !is.na(gene), .(
  genomewide_screened = TRUE,
  genomewide_cells = sum(num_cells_filtered, na.rm = TRUE),
  transcript_tests = sum(is.finite(energy_test_p_value)),
  genomewide_effect_p_min = min(energy_test_p_value, na.rm = TRUE),
  genomewide_fold_expr = weighted.mean(
    fold_expr, pmax(num_cells_filtered, 1), na.rm = TRUE
  )
), by = gene]
gw_gene[, genomewide_effect_p := pmin(
  genomewide_effect_p_min * pmax(transcript_tests, 1L), 1
)]
gw_gene[!is.finite(genomewide_effect_p), genomewide_effect_p := NA_real_]
gw_gene[, genomewide_effect_fdr := p.adjust(
  genomewide_effect_p, method = "BH"
)]

genes <- data.table(gene = sort(unique(c(
  rep_tasks$gene, norman_single$gene, combo_gene$gene, multi_gene$gene,
  gw_gene$gene
))))
evidence <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
  list(genes, rep_tasks, norman_single, combo_gene, multi_gene, gw_gene)
)
logical_cols <- c(
  "replogle_eligible", "norman_single_eligible", "multiperturb_eligible"
)
for (nm in logical_cols) evidence[is.na(get(nm)), (nm) := FALSE]
evidence[is.na(genomewide_screened), genomewide_screened := FALSE]
count_cols <- c(
  "observed_combinations", "synergistic_combinations",
  "antagonistic_combinations", "redundant_combinations",
  "norman_single_conditions", "multiperturb_contexts"
)
for (nm in count_cols) evidence[is.na(get(nm)), (nm) := 0]

evidence[, datasets_with_causal_support :=
  as.integer(replogle_eligible) +
  as.integer(norman_single_eligible) +
  as.integer(multiperturb_eligible)
]
evidence[, independent_contexts_with_causal_support :=
  as.integer(replogle_eligible | norman_single_eligible) +
  as.integer(multiperturb_eligible)
]
off_target_cut <- median(
  evidence$off_target_deviation[evidence$replogle_eligible],
  na.rm = TRUE
)
evidence[, path_specific :=
  (replogle_eligible & off_target_deviation <= off_target_cut) |
  (norman_single_eligible & norman_single_reproducibility >= 0.5)
]
evidence[, combination_consistent :=
  synergistic_combinations > antagonistic_combinations &
  observed_combinations >= 1L
]
evidence[, causal_regulator :=
  replogle_eligible | norman_single_eligible | multiperturb_eligible
]
evidence[, validated_control_node_candidate :=
  causal_regulator & independent_contexts_with_causal_support >= 2L &
  combination_consistent & path_specific
]
positive_combo_counts <- evidence[
  observed_combinations > 0L, observed_combinations
]
association_hub_cut <- if (length(positive_combo_counts)) {
  as.numeric(quantile(positive_combo_counts, 0.9, na.rm = TRUE))
} else {
  Inf
}
evidence[, evidence_class := fifelse(
  validated_control_node_candidate, "validated_control_node_candidate",
  fifelse(causal_regulator, "causal_regulator",
    fifelse(
      observed_combinations > 0L &
        observed_combinations >= association_hub_cut,
      "association_hub",
      fifelse(
        genomewide_screened & genomewide_effect_fdr < 0.05,
        "state_marker", "insufficient_causal_evidence"
      )
    )
  )
)]
evidence[, true_control_node := FALSE]
evidence[, true_control_node_reason :=
  "requires held-out cross-cell-line target reachability; not claimed here"
]
setorder(evidence, -validated_control_node_candidate,
         -datasets_with_causal_support, -synergistic_combinations,
         control_energy)
fwrite(evidence, file.path(out_dir, "control_node_evidence.tsv"), sep = "\t")

# Observed-support sequential policies. Only continuation of the same drug across
# measured 6 h and 24 h endpoints is evaluable; switches and ordering are abstained.
mix <- fread("results/mixseq_probabilities/pairs.tsv")
mix_time <- mix[time_hours %in% c(6, 24), .(
  p_apoptosis = weighted.mean(p_apoptosis, treated_cells),
  p_proliferation = weighted.mean(p_proliferation, treated_cells),
  p_stress = weighted.mean(p_stress, treated_cells),
  uncertainty = max(
    p_apoptosis_se, p_proliferation_se, p_stress_se, na.rm = TRUE
  ),
  treated_cells = sum(treated_cells)
), by = .(cell_line, perturbation, time_hours)]
wide <- dcast(mix_time, cell_line + perturbation ~ time_hours,
  value.var = c(
    "p_apoptosis", "p_proliferation", "p_stress",
    "uncertainty", "treated_cells"
  )
)
wide <- wide[complete.cases(wide[, .(
  p_apoptosis_6, p_apoptosis_24, p_proliferation_6,
  p_proliferation_24, p_stress_6, p_stress_24
)])]
wide[, early_role := fifelse(
  p_stress_6 >= pmax(p_apoptosis_6, p_proliferation_6),
  "early_stress_or_barrier_response", "early_state_movement"
)]
wide[, late_role := c("apoptosis", "proliferation", "stress")[
  max.col(as.matrix(.SD), ties.method = "first")
], .SDcols = c("p_apoptosis_24", "p_proliferation_24", "p_stress_24")]
wide[, observed_support_policy := fifelse(
  uncertainty_6 <= 0.10 & uncertainty_24 <= 0.10 &
    pmax(p_apoptosis_24, p_proliferation_24, p_stress_24) >
      pmax(p_apoptosis_6, p_proliferation_6, p_stress_6),
  "continue_same_intervention_to_24h",
  "abstain_no_supported_improvement"
)]
wide[, commitment_lock_status :=
  "abstain_no_washout_or_reversal_measurement"
]
fwrite(wide, file.path(out_dir, "observed_support_sequences.tsv"), sep = "\t")

sequence_boundary <- data.table(
  policy_space = c(
    "same_drug_6_to_24h", "switch_drug_A_to_B",
    "order_Norman_pair", "commitment_lock"
  ),
  status = c("evaluated", "abstain", "abstain", "abstain"),
  reason = c(
    "same context and intervention measured at both horizons",
    "no ordered switch observations",
    "pair is simultaneous and order is unobserved",
    "no washout or reversal endpoint"
  )
)
fwrite(sequence_boundary, file.path(out_dir, "sequence_support_boundary.tsv"),
       sep = "\t")

# Baseline execution registry prevents unavailable methods from being silently
# treated as completed comparisons.
baseline_status <- data.table(
  baseline = c(
    "differential_expression_ranking", "network_centrality",
    "SCENIC_regulon_ranking", "elastic_net", "GEARS", "CPA_chemCPA",
    "scGen", "TxPert", "static_shortest_path",
    "no_dynamics_classifier"
  ),
  status = c(
    "descriptive_ranking_only_not_a_matched_numeric_benchmark",
    "association_degree_only_not_network_centrality",
    "regulon_axes_used_as_covariates_not_a_ranking_benchmark",
    "implemented_in_MIXSeq_grouped_validation",
    "not_run_missing_compatible_checkpoint_and_training_runtime",
    "not_run_missing_compatible_checkpoint_and_training_runtime",
    "not_run_missing_compatible_checkpoint_and_training_runtime",
    "not_run_missing_compatible_checkpoint_and_training_runtime",
    "additive_response_null_not_a_static_shortest_path",
    "implemented_as_training_mean_and_ridge_static"
  )
)
baseline_status[, headline_comparison_allowed := baseline %in% c(
  "elastic_net", "no_dynamics_classifier"
)]
fwrite(baseline_status, file.path(out_dir, "baseline_status.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "target_basins", "fragile_target_basins", "redundant_target_basins",
    "causal_regulators", "validated_control_node_candidates",
    "true_control_nodes_claimed", "supported_same_drug_sequences"
  ),
  value = c(
    nrow(degeneracy),
    sum(degeneracy$fragile_control_axis),
    sum(degeneracy$redundant_control_landscape),
    sum(evidence$causal_regulator),
    sum(evidence$validated_control_node_candidate),
    sum(evidence$true_control_node),
    sum(wide$observed_support_policy ==
          "continue_same_intervention_to_24h")
  )
)
fwrite(summary, file.path(out_dir, "summary.tsv"), sep = "\t")
print(summary)
print(head(evidence[validated_control_node_candidate == TRUE], 20))
