#!/usr/bin/env Rscript

# Task-registry contract
# A task is a complete observed perturbation condition with a matched source
# distribution and a reproducible endpoint distribution. Endpoint-response basins
# and phenotype-labelled fate tasks are kept as separate evidence tiers.

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
})
source("code/common.R")

set.seed(1903)
out_dir <- "results/task_registry"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

program_cols <- c(
  apoptosis = "score_apoptosis",
  proliferation = "score_proliferation",
  stress = "score_stress"
)

# Replogle K562: program-defined endpoints with independent gem-group support.
rep_cells <- as.data.table(readRDS("inputs/replogle_scores.rds"))
rep_cells[, gene := as.character(gene)]
rep_cells <- rep_cells[
  passes_qc_final == TRUE & scDblFinder.class == "singlet"
]
rep_cells[, `:=`(
  log_umi = log1p(UMI_count),
  log_genes = log1p(n_genes_by_counts)
)]
rep_control <- rep_cells[gene == "non-targeting"]
control_mean <- vapply(program_cols, function(v) {
  mean(rep_control[[v]], na.rm = TRUE)
}, numeric(1))
control_sd <- vapply(program_cols, function(v) {
  sd(rep_control[[v]], na.rm = TRUE)
}, numeric(1))
control_q10 <- vapply(program_cols, function(v) {
  quantile(rep_control[[v]], 0.10, na.rm = TRUE)
}, numeric(1))
control_q90 <- vapply(program_cols, function(v) {
  quantile(rep_control[[v]], 0.90, na.rm = TRUE)
}, numeric(1))

rep_summary <- rep_cells[gene != "non-targeting", {
  means <- vapply(program_cols, function(v) mean(get(v), na.rm = TRUE), numeric(1))
  delta_z <- (means - control_mean) / pmax(control_sd, 1e-8)
  primary_idx <- which.max(abs(delta_z))
  primary_program <- names(program_cols)[primary_idx]
  primary_direction <- ifelse(delta_z[primary_idx] >= 0, "up", "down")
  primary_values <- get(program_cols[[primary_program]])
  nuisance_means <- c(
    log_umi = mean(log_umi, na.rm = TRUE),
    mitopercent = mean(mitopercent, na.rm = TRUE),
    log_genes = mean(log_genes, na.rm = TRUE)
  )
  success <- if (primary_direction == "up") {
    mean(primary_values > control_q90[[primary_program]], na.rm = TRUE)
  } else {
    mean(primary_values < control_q10[[primary_program]], na.rm = TRUE)
  }
  list(
    target_cell_count = .N,
    n_gem_groups = uniqueN(gem_group),
    apoptosis_delta_z = delta_z[["apoptosis"]],
    proliferation_delta_z = delta_z[["proliferation"]],
    stress_delta_z = delta_z[["stress"]],
    target_state_id = paste0(
      "K562_", primary_program, "_", primary_direction
    ),
    primary_program = primary_program,
    primary_direction = primary_direction,
    target_state_distance = safe_norm(delta_z),
    control_success_rate = success,
    failure_probability = 1 - success,
    off_target_deviation = safe_norm(delta_z[-primary_idx]),
    control_energy = 1 / (sum(delta_z^2) + 1e-8),
    nuisance_log_umi = nuisance_means[["log_umi"]],
    nuisance_mitopercent = nuisance_means[["mitopercent"]],
    nuisance_log_genes = nuisance_means[["log_genes"]]
  )
}, by = gene]
control_nuisance_mean <- c(
  log_umi = mean(rep_control$log_umi, na.rm = TRUE),
  mitopercent = mean(rep_control$mitopercent, na.rm = TRUE),
  log_genes = mean(rep_control$log_genes, na.rm = TRUE)
)
control_nuisance_sd <- c(
  log_umi = sd(rep_control$log_umi, na.rm = TRUE),
  mitopercent = sd(rep_control$mitopercent, na.rm = TRUE),
  log_genes = sd(rep_control$log_genes, na.rm = TRUE)
)
control_nuisance_sd[!is.finite(control_nuisance_sd) |
  control_nuisance_sd < 1e-8] <- 1
rep_summary[, nuisance_max_abs_z := pmax(
  abs((nuisance_log_umi - control_nuisance_mean[["log_umi"]]) /
    control_nuisance_sd[["log_umi"]]),
  abs((nuisance_mitopercent - control_nuisance_mean[["mitopercent"]]) /
    control_nuisance_sd[["mitopercent"]]),
  abs((nuisance_log_genes - control_nuisance_mean[["log_genes"]]) /
    control_nuisance_sd[["log_genes"]])
)]
rep_summary[, nuisance_separation_check :=
  is.finite(nuisance_max_abs_z) &
  target_state_distance > nuisance_max_abs_z
]

rep_by_gem <- rep_cells[, lapply(.SD, mean, na.rm = TRUE),
  by = .(gene, gem_group), .SDcols = unname(program_cols)
]
control_by_gem <- rep_by_gem[gene == "non-targeting"]
control_by_gem[, gene := NULL]
setnames(control_by_gem, unname(program_cols), paste0(names(program_cols), "_ctrl"))
rep_by_gem <- merge(
  rep_by_gem[gene != "non-targeting"], control_by_gem,
  by = "gem_group", all.x = TRUE
)
for (nm in names(program_cols)) {
  rep_by_gem[, paste0(nm, "_dz") :=
    (get(program_cols[[nm]]) - get(paste0(nm, "_ctrl"))) /
      pmax(control_sd[[nm]], 1e-8)]
}

rep_repro <- rep_by_gem[, {
  overall <- as.numeric(rep_summary[match(.BY$gene, gene),
    .(apoptosis_delta_z, proliferation_delta_z, stress_delta_z)
  ])
  mat <- as.matrix(.SD)
  cosines <- apply(mat, 1, safe_cosine, y = overall)
  selection_label <- vapply(seq_len(nrow(mat)), function(i) {
    train_mean <- colMeans(mat[-i, , drop = FALSE], na.rm = TRUE)
    primary <- which.max(abs(train_mean))
    paste0(
      names(program_cols)[primary], "_",
      ifelse(train_mean[primary] >= 0, "up", "down")
    )
  }, character(1))
  held_signed_effect <- vapply(seq_len(nrow(mat)), function(i) {
    train_mean <- colMeans(mat[-i, , drop = FALSE], na.rm = TRUE)
    primary <- which.max(abs(train_mean))
    direction <- ifelse(train_mean[primary] >= 0, 1, -1)
    mat[i, primary] * direction
  }, numeric(1))
  valid_effect <- held_signed_effect[is.finite(held_signed_effect)]
  selected_target <- names(which.max(table(selection_label)))
  se <- standard_error(valid_effect)
  df <- length(valid_effect) - 1L
  effect_mean <- mean(valid_effect)
  statistic <- effect_mean / se
  p_value <- if (df >= 1L && is.finite(statistic)) {
    2 * pt(-abs(statistic), df = df)
  } else if (df >= 1L && effect_mean > 0 && se == 0) {
    0
  } else {
    NA_real_
  }
  critical <- if (df >= 1L) qt(0.975, df = df) else NA_real_
  list(
    replicate_groups = nrow(mat),
    median_replicate_cosine = median(cosines, na.rm = TRUE),
    selected_target_state = selected_target,
    target_selection_consistency = max(table(selection_label)) /
      length(selection_label),
    concordant_replicates = sum(valid_effect > 0),
    primary_effect_mean = effect_mean,
    primary_effect_se = se,
    primary_effect_df = df,
    primary_effect_p = p_value,
    primary_ci_low = effect_mean - critical * se,
    primary_ci_high = effect_mean + critical * se,
    control_success_rate_crossfit = mean(valid_effect > 0),
    target_selection_method =
      "leave_one_gem_group_out_program_and_direction"
  )
}, by = gene, .SDcols = paste0(names(program_cols), "_dz")]
rep_summary <- merge(rep_summary, rep_repro, by = "gene", all.x = TRUE)
rep_summary[, primary_effect_fdr := p.adjust(primary_effect_p, method = "BH")]
rep_summary[, `:=`(
  primary_program = sub("_(up|down)$", "", selected_target_state),
  primary_direction = sub("^.*_(up|down)$", "\\1", selected_target_state),
  target_state_id = paste0("K562_", selected_target_state),
  control_success_rate = control_success_rate_crossfit,
  failure_probability = 1 - control_success_rate_crossfit,
  uncertainty = primary_effect_se
)]
rep_summary[, off_target_deviation := mapply(
  function(a, p, s, primary) {
    values <- c(apoptosis = a, proliferation = p, stress = s)
    safe_norm(values[names(values) != primary])
  },
  apoptosis_delta_z, proliferation_delta_z, stress_delta_z,
  primary_program
)]
rep_summary[, passes_task_qc :=
  target_cell_count >= 100L &
  replicate_groups >= 3L &
  target_selection_consistency >= 0.6 &
  concordant_replicates >= ceiling(0.6 * replicate_groups) &
  primary_ci_low > 0 &
  is.finite(primary_effect_fdr) &
  primary_effect_fdr < 0.05 &
  nuisance_separation_check
]
fwrite(
  rep_summary[, .(
    gene, target_cell_count, replicate_groups, selected_target_state,
    target_selection_consistency, concordant_replicates,
    primary_effect_mean, primary_effect_se, primary_effect_df,
    primary_effect_p, primary_effect_fdr, primary_ci_low, primary_ci_high,
    control_success_rate_crossfit, median_replicate_cosine,
    nuisance_max_abs_z, nuisance_separation_check, passes_task_qc,
    target_selection_method
  )],
  file.path(out_dir, "replogle_task_qc.tsv"),
  sep = "\t"
)

rep_tasks <- rep_summary[, .(
  task_id = paste0("replogle_k562|", gene),
  dataset_id = "replogle_2022_k562_essential",
  task_type = "phenotype_program_endpoint",
  evidence_tier = "single_cell_perturbation",
  source_state_id = "K562_non_targeting",
  target_state_id,
  context_id = "K562_CRISPRi",
  horizon = NA_real_,
  perturbation_identity = gene,
  perturbation_type = "genetic_CRISPRi",
  intervention_size = 1L,
  source_cell_count = nrow(rep_control),
  target_cell_count,
  target_definition_source =
    "cross_fitted_selection_within_prespecified_program_set",
  target_reproducibility = pmin(
    median_replicate_cosine, target_selection_consistency
  ),
  nuisance_separation_check,
  target_state_distance,
  control_success_rate,
  failure_probability,
  off_target_deviation,
  control_energy,
  uncertainty,
  status = fifelse(passes_task_qc, "eligible", "excluded"),
  exclusion_reason = fifelse(passes_task_qc, NA_character_,
    "cell_count_replicate_or_effect_CI_gate"
  )
)]

# Norman: endpoint-response basins learned from full normalized response vectors.
norman_response <- as.data.table(readRDS("results/norman_control/condition_response_delta.rds"))
norman_index <- fread("results/norman_control/condition_index.tsv")
norman_qc <- fread("results/norman_control/condition_qc.tsv")
norman <- merge(norman_index, norman_qc,
  by = c("condition_id", "raw_guide_identity", "n_cells", "condition_type")
)
feature_cols <- setdiff(names(norman_response), "condition_id")
mat <- as.matrix(norman_response[, ..feature_cols])
row.names(mat) <- norman_response$condition_id
feature_sd <- apply(mat, 2, sd, na.rm = TRUE)
keep_feature <- is.finite(feature_sd) & feature_sd > 1e-6
pc <- prcomp(mat[, keep_feature, drop = FALSE], rank. = 20, scale. = TRUE)
perturbed_idx <- which(norman$condition_type != "control")
pc_perturbed <- pc$x[match(
  norman$condition_id[perturbed_idx], rownames(pc$x)
), , drop = FALSE]
k_candidates <- 8:24
silhouette_mean <- vapply(k_candidates, function(k) {
  km <- kmeans(pc_perturbed, centers = k, nstart = 20, iter.max = 100)
  mean(cluster::silhouette(km$cluster, dist(pc_perturbed))[, "sil_width"])
}, numeric(1))
best_k <- k_candidates[which.max(silhouette_mean)]
km <- kmeans(pc_perturbed, centers = best_k, nstart = 100, iter.max = 200)
basin_map <- data.table(
  condition_id = norman$condition_id[perturbed_idx],
  target_state_id = sprintf("NORMAN_BASIN_%02d", km$cluster)
)
norman <- merge(norman, basin_map, by = "condition_id", all.x = TRUE)
norman_source_n <- sum(norman$n_cells[norman$condition_type == "control"])
norman[, passes_registry_qc :=
  condition_type != "control" & n_cells >= 100L & passes_task_qc == TRUE
]
norman_tasks <- norman[condition_type != "control", .(
  task_id = paste0("norman_k562|", condition_id),
  dataset_id = "norman_2019_combinatorial",
  task_type = "observed_response_basin",
  evidence_tier = "single_cell_perturbation",
  source_state_id = "K562_negative_control",
  target_state_id,
  context_id = "K562_CRISPRa",
  horizon = NA_real_,
  perturbation_identity = condition_id,
  perturbation_type = "genetic_CRISPRa",
  intervention_size = fifelse(condition_type == "combination", 2L, 1L),
  source_cell_count = norman_source_n,
  target_cell_count = n_cells,
  target_definition_source = paste0(
    "unsupervised_response_basin_k", best_k,
    "_not_a_phenotype_label"
  ),
  target_reproducibility = median_replicate_cosine,
  nuisance_separation_check,
  target_state_distance = response_norm,
  control_success_rate = NA_real_,
  failure_probability = NA_real_,
  off_target_deviation = NA_real_,
  control_energy = fifelse(condition_type == "combination", 2, 1) /
    (response_norm^2 + 1e-8),
  uncertainty = 1 - pmax(pmin(median_replicate_cosine, 1), -1),
  status = fifelse(passes_registry_qc, "eligible", "excluded"),
  exclusion_reason = fifelse(passes_registry_qc, NA_character_,
    "minimum_100_cells_or_replicate_concordance"
  )
)]

fwrite(data.table(k = k_candidates, mean_silhouette = silhouette_mean),
  file.path(out_dir, "norman_basin_selection.tsv"), sep = "\t"
)

# MIX-Seq: matched DMSO-to-drug program endpoints.
mix <- fread("results/mixseq_probabilities/pairs.tsv")
mix[, primary_program := c("apoptosis", "proliferation", "stress")[
  max.col(as.matrix(.SD), ties.method = "first")
], .SDcols = c("p_apoptosis", "p_proliferation", "p_stress")]
mix[, target_state_id := paste0("drug_induced_", primary_program)]
mix[, target_state_distance := sqrt(
  (p_apoptosis - 0.1)^2 + (p_proliferation - 0.1)^2 +
    (p_stress - 0.1)^2
)]
mix[, control_success_rate := pmax(p_apoptosis, p_proliferation, p_stress)]
mix[, uncertainty := pmax(
  p_apoptosis_se, p_proliferation_se, p_stress_se, na.rm = TRUE
)]
mix[, passes_task_qc :=
  baseline_cells >= 100L & treated_cells >= 100L &
  is.finite(uncertainty) & uncertainty <= 0.10
]
mix_tasks <- mix[, .(
  task_id = pair_id,
  dataset_id = "mixseq_2020",
  task_type = "phenotype_program_endpoint",
  evidence_tier = "derived_matched_control",
  source_state_id = paste0(
    cell_line, "_", baseline, "_experiment_", experiment
  ),
  target_state_id,
  context_id = paste0(cell_line, "_experiment_", experiment),
  horizon = time_hours,
  perturbation_identity = perturbation,
  perturbation_type = task,
  intervention_size = 1L,
  source_cell_count = baseline_cells,
  target_cell_count = treated_cells,
  target_definition_source = "matched_control_program_tail_probability",
  target_reproducibility = 1 - uncertainty,
  nuisance_separation_check = NA,
  target_state_distance,
  control_success_rate,
  failure_probability = 1 - control_success_rate,
  off_target_deviation = rowSums(as.matrix(.SD)) - control_success_rate,
  control_energy = 1 / (target_state_distance^2 + 1e-8),
  uncertainty,
  status = fifelse(passes_task_qc, "eligible", "excluded"),
  exclusion_reason = fifelse(passes_task_qc, NA_character_,
    "minimum_100_cells_or_probability_SE"
  )
), .SDcols = c("p_apoptosis", "p_proliferation", "p_stress")]

# MultiPerturb-seq: multi-guide reproducible joint-modality endpoints.
multi <- fread("results/multiperturb_multiome/target_reproducibility.tsv")
multi_tasks <- multi[, .(
  task_id = paste0("multiperturb|", context, "|", guide_gene),
  dataset_id = "yan_2024_multiperturbseq",
  task_type = "joint_modality_endpoint",
  evidence_tier = "single_cell_perturb_multiome",
  source_state_id = paste0(context, "_non_targeting"),
  target_state_id = paste0(context, "_observed_joint_endpoint"),
  context_id = context,
  horizon = NA_real_,
  perturbation_identity = guide_gene,
  perturbation_type = "genetic_CRISPR",
  intervention_size = 1L,
  source_cell_count = source_control_cells,
  target_cell_count = n_cells,
  target_definition_source = "joint_RNA_ATAC_response",
  target_reproducibility = pmin(
    median_rna_guide_cosine, median_atac_guide_cosine
  ),
  nuisance_separation_check,
  target_state_distance = sqrt(rna_movement^2 + atac_movement^2),
  control_success_rate = NA_real_,
  failure_probability = NA_real_,
  off_target_deviation = abs(rna_movement - atac_movement),
  control_energy = 1 / (rna_movement^2 + atac_movement^2 + 1e-8),
  uncertainty = 1 - pmin(
    median_rna_guide_cosine, median_atac_guide_cosine
  ),
  status = fifelse(
    passes_task_qc & n_cells >= 100L & source_control_cells >= 100L,
    "eligible", "excluded"
  ),
  exclusion_reason = fifelse(
    passes_task_qc & n_cells >= 100L & source_control_cells >= 100L,
    NA_character_,
    "minimum_cells_multi_guide_reproducibility_or_nuisance_gate"
  )
)]

registry <- rbindlist(list(rep_tasks, norman_tasks, mix_tasks, multi_tasks),
  use.names = TRUE, fill = TRUE
)
registry[, modality_mask := c(
  replogle_2022_k562_essential =
    "RNA|pathway_programs|cell_cycle|metabolic_proxy",
  norman_2019_combinatorial = "RNA",
  mixseq_2020 = "RNA|baseline_latent|pathway_programs|cell_cycle",
  yan_2024_multiperturbseq = "RNA|ATAC"
)[dataset_id]]
registry[, control_energy_definition :=
  "relative_effort_proxy_intervention_size_over_squared_response_distance"
]
registry[, atomic_split_unit := paste(
  dataset_id, perturbation_identity, context_id,
  fifelse(is.na(horizon), "endpoint", as.character(horizon)),
  sep = "|"
)]
if (anyDuplicated(registry$atomic_split_unit)) {
  duplicate_units <- registry[
    duplicated(atomic_split_unit) |
      duplicated(atomic_split_unit, fromLast = TRUE),
    unique(atomic_split_unit)
  ]
  stop(
    "Non-unique atomic task units: ",
    paste(duplicate_units, collapse = ", ")
  )
}
setorder(registry, dataset_id, task_id)
fwrite(registry, file.path(out_dir, "state_transition_tasks.tsv"), sep = "\t")

gate_summary <- registry[, .(
  candidate_tasks = .N,
  eligible_tasks = sum(status == "eligible"),
  phenotype_labelled_tasks = sum(
    status == "eligible" & task_type == "phenotype_program_endpoint"
  ),
  observed_response_basin_tasks = sum(
    status == "eligible" & task_type == "observed_response_basin"
  ),
  median_target_cells = as.numeric(median(target_cell_count, na.rm = TRUE))
), by = dataset_id]
fwrite(gate_summary, file.path(out_dir, "task_gate_by_dataset.tsv"), sep = "\t")

eligible_n<-registry[status=="eligible",.N]
cat("Norman response basins:", best_k, "\n")
print(gate_summary)
cat("\nTotal eligible tasks:", eligible_n, "\n")
