#!/usr/bin/env Rscript
source("code/common.R")
out_dir<-"results/norman_control";dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
x<-read_norman_input();aux<-readRDS("inputs/auxiliary.rds")
delta<-x$delta;delta_gem<-x$delta_gem;feature_names<-x$features;selected<-seq_along(feature_names)
guide_simple<-x$conditions;n_condition<-length(guide_simple);n_gem<-dim(delta_gem)[3];n_condition_qc<-x$n_cells
n_condition_gem<-as.vector(x$n_condition_gem);is_single<-x$is_single;is_pair<-x$is_pair;is_control<-x$is_control;single_gene<-x$single_gene
raw_identity_by_condition<-aux$norman_index$raw_guide_identity[match(guide_simple,aux$norman_index$condition_id)]
total_condition<-aux$norman_qc$mean_total_counts[match(guide_simple,aux$norman_qc$condition_id)]*n_condition_qc
single_indices <- split(which(is_single & !is.na(single_gene)), single_gene[
  is_single & !is.na(single_gene)
])
single_delta <- lapply(single_indices, function(idx) {
  w <- n_condition_qc[idx]
  as.numeric(delta[, idx, drop = FALSE] %*% w / sum(w))
})
single_delta_gem <- lapply(single_indices, function(idx) {
  out <- matrix(NA_real_, nrow = length(selected), ncol = n_gem)
  for (g in seq_len(n_gem)) {
    w <- n_condition_gem[idx + (g - 1L) * n_condition]
    if (sum(w) > 0) {
      slice <- matrix(
        delta_gem[, idx, g, drop = FALSE],
        nrow = length(selected), ncol = length(idx)
      )
      out[, g] <- as.numeric(slice %*% w / sum(w))
    }
  }
  out
})

condition_qc <- rbindlist(lapply(seq_len(n_condition), function(j) {
  d <- delta[, j]
  gem_cos <- vapply(seq_len(n_gem), function(g) {
    safe_cosine(delta_gem[, j, g], d)
  }, numeric(1))
  data.table(
    condition_id = guide_simple[j],
    raw_guide_identity = raw_identity_by_condition[j],
    condition_type = fifelse(is_pair[j], "combination",
      fifelse(is_single[j], "single", "control")
    ),
    n_cells = n_condition_qc[j],
    n_gem_groups = sum(n_condition_gem[
      j + (seq_len(n_gem) - 1L) * n_condition
    ] > 0L),
    response_norm = safe_norm(d),
    median_replicate_cosine = median(gem_cos, na.rm = TRUE),
    concordant_replicates = sum(gem_cos > 0.2, na.rm = TRUE),
    mean_total_counts = total_condition[j] / pmax(n_condition_qc[j], 1L),
    passes_task_qc = n_condition_qc[j] >= 30L &
      sum(gem_cos > 0.2, na.rm = TRUE) >= 6L
  )
}))
control_log_counts <- log1p(
  condition_qc$mean_total_counts[condition_qc$condition_type == "control"]
)
control_log_counts_sd <- sd(control_log_counts, na.rm = TRUE)
if (!is.finite(control_log_counts_sd) || control_log_counts_sd < 1e-8) {
  control_log_counts_sd <- 1
}
condition_qc[, library_size_z := (
  log1p(mean_total_counts) - mean(control_log_counts, na.rm = TRUE)
) / control_log_counts_sd]
condition_qc[, nuisance_separation_check :=
  is.finite(library_size_z) & response_norm > abs(library_size_z)
]
condition_qc[, passes_task_qc := passes_task_qc &
  nuisance_separation_check
]

pair_rows <- which(is_pair)
interaction <- rbindlist(lapply(pair_rows, function(j) {
  genes <- strsplit(guide_simple[j], "_", fixed = TRUE)[[1]]
  if (length(genes) != 2L || !all(genes %in% names(single_delta))) return(NULL)
  da <- single_delta[[genes[1]]]
  db <- single_delta[[genes[2]]]
  dp <- delta[, j]
  additive <- da + db
  int <- dp - additive
  # The tested direction must not be defined by the pair response itself.
  # Singles define the expected target direction; the pair is then tested
  # against that independent additive direction.
  target_unit <- additive / pmax(safe_norm(additive), 1e-10)
  projection <- sum(int * target_unit)
  orthogonal <- int - projection * target_unit

  projection_g <- vapply(seq_len(n_gem), function(g) {
    ig <- delta_gem[, j, g] -
      single_delta_gem[[genes[1]]][, g] -
      single_delta_gem[[genes[2]]][, g]
    sum(ig * target_unit, na.rm = TRUE)
  }, numeric(1))
  valid_projection <- projection_g[is.finite(projection_g)]
  n_projection <- length(valid_projection)
  projection_mean <- mean(valid_projection)
  se <- standard_error(valid_projection)
  df <- n_projection - 1L
  t_statistic <- projection_mean / se
  p_value <- if (
    n_projection >= 3L && is.finite(t_statistic)
  ) {
    2 * pt(-abs(t_statistic), df = df)
  } else if (
    n_projection >= 3L && is.finite(projection_mean) &&
      projection_mean != 0 && isTRUE(all.equal(se, 0))
  ) {
    0
  } else {
    NA_real_
  }
  critical <- if (df >= 1L) qt(0.975, df = df) else NA_real_
  ci_low <- projection_mean - critical * se
  ci_high <- projection_mean + critical * se
  movement_a <- safe_norm(da)
  movement_b <- safe_norm(db)
  movement_pair <- safe_norm(dp)
  redundant <- max(safe_cosine(dp, da), safe_cosine(dp, db), na.rm = TRUE) > 0.8 &&
    movement_pair <= 1.1 * max(movement_a, movement_b)
  data.table(
    combination = guide_simple[j],
    gene_a = genes[1],
    gene_b = genes[2],
    n_cells = n_condition_qc[j],
    target_movement = movement_pair,
    single_a_movement = movement_a,
    single_b_movement = movement_b,
    additive_movement = safe_norm(additive),
    interaction_norm = safe_norm(int),
    interaction_projection = projection,
    replicate_projection_mean = projection_mean,
    replicate_projection_se = se,
    replicate_projection_df = df,
    replicate_projection_t = t_statistic,
    replicate_projection_p = p_value,
    ci_low = ci_low,
    ci_high = ci_high,
    additive_cosine = safe_cosine(dp, additive),
    off_target_interaction = safe_norm(orthogonal),
    redundant_candidate = redundant,
    control_energy = 2 / (movement_pair^2 + 1e-8)
  )
}))
interaction[, replicate_projection_fdr := p.adjust(
  replicate_projection_p, method = "BH"
)]
interaction[, interaction_class := fifelse(
  is.finite(replicate_projection_fdr) &
    replicate_projection_fdr < 0.05 &
    replicate_projection_mean > 0,
  "synergistic",
  fifelse(
    is.finite(replicate_projection_fdr) &
      replicate_projection_fdr < 0.05 &
      replicate_projection_mean < 0,
    "antagonistic",
    fifelse(redundant_candidate, "redundant", "context-dependent")
  )
)]
interaction[, minimum_intervention_size := fifelse(
  interaction_class == "synergistic" &
    target_movement > 1.1 * pmax(single_a_movement, single_b_movement),
  2L, 1L
)]

condition_index <- data.table(
  condition_id = guide_simple,
  raw_guide_identity = raw_identity_by_condition,
  n_cells = n_condition_qc,
  condition_type = fifelse(is_pair, "combination",
    fifelse(is_single, "single", "control")
  ),
  single_gene = single_gene
)

response_export <- as.data.table(t(delta))
setnames(response_export, feature_names)
response_export[, condition_id := guide_simple]
setcolorder(response_export, "condition_id")

fwrite(condition_index, file.path(out_dir, "condition_index.tsv"), sep = "\t")
fwrite(condition_qc, file.path(out_dir, "condition_qc.tsv"), sep = "\t")
fwrite(interaction, file.path(out_dir, "combination_interactions.tsv"), sep = "\t")
saveRDS(response_export,file.path(out_dir,"condition_response_delta.rds"),compress="xz")

summary <- rbindlist(list(
  data.table(metric = "cells_after_metadata_qc", value = sum(n_condition_qc)),
  data.table(metric = "observed_conditions", value = n_condition),
  data.table(metric = "single_conditions", value = sum(is_single)),
  data.table(metric = "combination_conditions", value = sum(is_pair)),
  data.table(metric = "eligible_conditions", value = sum(condition_qc$passes_task_qc)),
  data.table(metric = "synergistic_pairs",
             value = sum(interaction$interaction_class == "synergistic")),
  data.table(metric = "antagonistic_pairs",
             value = sum(interaction$interaction_class == "antagonistic")),
  data.table(metric = "redundant_pairs",
             value = sum(interaction$interaction_class == "redundant")),
  data.table(metric = "context_dependent_pairs",
             value = sum(interaction$interaction_class == "context-dependent"))
))
fwrite(summary, file.path(out_dir, "summary.tsv"), sep = "\t")
print(summary)
