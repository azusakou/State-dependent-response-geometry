#!/usr/bin/env Rscript

# Main-figure contracts
# Fig. 2: baseline state predicts unequal local reachable fate dimensions.
# Fig. 3: model performance depends on the held-out axis and abstention support.
# Fig. 4: observed combinations show reproducible non-additive control.
# Fig. 5: perturb-multiome separates transcription- and chromatin-dominant movers.
# Fig. 6: target basins differ in control degeneracy and substitutability.
# Fig. 7: only same-drug continuation is supported for virtual sequences.
# Fig. 8: the atlas distinguishes causal evidence from true-control-node claims.
# Archetypes and panel-specific review risks: docs/figure_redesign_contract.md.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("R/figure_theme.R")
theme_set(theme_nature_control())

source_dir <- "results/figure_source_data"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures/main", recursive = TRUE, showWarnings = FALSE)

method_colors <- c(
  training_mean = "#A0A0A0",
  ridge_static = "#7A7A7A",
  elastic_net = "#4C78A8",
  local_response_ridge = "#E3A857",
  kernel_transition = "#2A9D8F",
  conditional_rate_network = "#7B6BA8",
  kernel_regression = "#2A9D8F",
  additive_null = "#A0A0A0",
  nonlinear_xgboost = "#C95A4B"
)
class_colors <- c(
  synergistic = "#2A9D8F",
  antagonistic = "#C95A4B",
  redundant = "#7B6BA8",
  `context-dependent` = "#E3A857"
)

# Figure 2 -----------------------------------------------------------------
land_unit <- fread("results/larry_controllability/local_controllability.tsv")
if ("local_effective_rank_stress_state_loo" %in% names(land_unit)) {
  land_unit[, local_effective_rank :=
    local_effective_rank_stress_state_loo
  ]
}
land_unit[, dataset_label := fifelse(
  dataset_id == "larry_cytokine", "Cytokine perturbation", "In-vitro lineage"
)]
state_features <- c(
  "cell_cycle_state", "stress_state", "chromatin_regulator_proxy",
  "metabolic_state", "mitochondrial_fraction"
)
# A clone can appear under more than one future condition. Figures use one
# equally weighted row per biological clone; the condition-specific rows remain
# available in the source table used by the clustered regression.
land <- land_unit[, lapply(.SD, mean, na.rm = TRUE),
  by = .(dataset_id, dataset_label, clone_id),
  .SDcols = c(
    state_features, "local_effective_rank", "local_transition_entropy"
  )
]
land[, stress_bin := cut(
  stress_state,
  breaks = unique(quantile(stress_state, seq(0, 1, length.out = 21),
                           na.rm = TRUE)),
  include.lowest = TRUE, labels = FALSE
), by = dataset_id]
stress_summary <- land[, .(
  stress_state = mean(stress_state, na.rm = TRUE),
  effective_rank = mean(local_effective_rank, na.rm = TRUE),
  se = sd(local_effective_rank, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(dataset_id, dataset_label, stress_bin)]
land[, stress_quartile := factor(
  pmin(4L, pmax(1L, ceiling(4 * frank(
    stress_state, ties.method = "average"
  ) / .N))),
  levels = 1:4, labels = c("Q1 low", "Q2", "Q3", "Q4 high")
), by = dataset_id]
quartile_summary <- land[, .(
  effective_rank = mean(local_effective_rank, na.rm = TRUE),
  se = sd(local_effective_rank, na.rm = TRUE) / sqrt(.N)
), by = .(dataset_label, stress_quartile)]
assoc <- fread("results/larry_controllability/state_associations.tsv")[
  outcome == "local_effective_rank" & term %in% state_features
]
assoc[, signed_evidence := sign(estimate) * -log10(pmax(p_adj_bh, 1e-50))]
feature_sd <- melt(land[, c("dataset_id", state_features), with = FALSE],
  id.vars = "dataset_id", variable.name = "term")[,
    .(predictor_sd = sd(value, na.rm = TRUE)), by = .(dataset_id, term)]
assoc <- merge(assoc, feature_sd, by = c("dataset_id", "term"))
assoc[, `:=`(scaled_estimate = estimate * predictor_sd,
              scaled_se = std_error * predictor_sd)]
assoc[, term_label := c(
  cell_cycle_state = "Cell cycle",
  stress_state = "Stress",
  chromatin_regulator_proxy = "Chromatin-regulator proxy",
  metabolic_state = "Oxidative metabolism",
  mitochondrial_fraction = "Mitochondrial fraction"
)[term]]

p2a <- ggplot(stress_summary, aes(stress_state, effective_rank,
                                  colour = dataset_label)) +
  geom_errorbar(aes(ymin = effective_rank - se, ymax = effective_rank + se),
                width = 0, linewidth = 0.3, alpha = 0.6) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c(
    "Cytokine perturbation" = palette_control[["genetic"]],
    "In-vitro lineage" = palette_control[["dynamics"]]
  )) +
  labs(
    title = "Stress associates with local fate-response dimensionality",
    x = "Baseline stress-program activity",
    y = "Local response effective rank",
    colour = NULL
  ) +
  theme(legend.position = "top")

p2b <- ggplot(land, aes(local_transition_entropy, local_effective_rank)) +
  geom_bin_2d(bins = 18) +
  facet_wrap(~dataset_label, strip.position = "bottom") +
  scale_x_continuous(breaks = c(0, 0.75, 1.5)) +
  scale_fill_gradient(low = "#EAF2F6", high = "#356C91", trans = "sqrt") +
  labs(x = "Local transition entropy", y = "Effective rank", fill = "Clones") +
  theme(legend.position = "bottom")

p2c <- ggplot(assoc, aes(scaled_estimate, term_label,
                         colour = dataset_id)) +
  geom_vline(xintercept = 0, colour = palette_control[["neutral_light"]]) +
  geom_errorbar(aes(xmin = scaled_estimate - 1.96 * scaled_se,
                    xmax = scaled_estimate + 1.96 * scaled_se),
                orientation = "y", width = 0, linewidth = 0.4,
                position = position_dodge(width = 0.45)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c(
    larry_cytokine = palette_control[["genetic"]],
    larry_in_vitro = palette_control[["dynamics"]]
  ), labels = c("Cytokine perturbation", "In-vitro lineage")) +
  labs(
    title = "State covariates associate with local response geometry",
    subtitle = "Signed evidence: sign(effect) × -log10(BH-adjusted P)",
    x = "Rank change per predictor SD (approximate 95% CI)",
    y = NULL, colour = NULL
  ) +
  theme(legend.position = "none")

p2d <- ggplot(land, aes(local_effective_rank, fill = dataset_label)) +
  geom_density(alpha = 0.35, linewidth = 0.45) +
  scale_fill_manual(values = c(
    "Cytokine perturbation" = palette_control[["genetic"]],
    "In-vitro lineage" = palette_control[["dynamics"]]
  )) +
  labs(
    title = "Response dimensionality varies among baseline clones",
    x = "Local response effective rank",
    y = "Density", fill = NULL
  ) +
  theme(legend.position = "none")

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(stress_summary, file.path(source_dir, "Fig2a_stress_bins.tsv"), sep = "\t")
fwrite(quartile_summary, file.path(source_dir, "Fig2b_stress_quartiles.tsv"),
       sep = "\t")
fwrite(assoc, file.path(source_dir, "Fig2c_state_associations.tsv"), sep = "\t")
fwrite(land, file.path(source_dir, "Fig2bd_clone_geometry.tsv"), sep = "\t")
save_pub_r(fig2, "figures/main/Fig2_local_controllability",
           width_mm = 183, height_mm = 145)

# Figure 3 -----------------------------------------------------------------
model <- fread("results/mixseq_control_models/metrics.tsv")
model <- model[split %in% c(
  "leave_cell_line", "leave_target_state", "leave_perturbation"
)]
model[, split_label := c(
  leave_cell_line = "Cell line",
  leave_target_state = "Target state",
  leave_perturbation = "Perturbation"
)[split]]
model[, method_label := gsub("_", " ", method)]
model[, method_label := c(
  training_mean = "training mean",
  ridge_static = "ridge static",
  elastic_net = "elastic net",
  local_response_ridge = "local ridge",
  kernel_transition = "kernel transition",
  conditional_rate_network = "conditional rate"
)[method]]

model[, displayed_rmse := fifelse(coverage > 0, target_state_rmse, NA_real_)]
p3a <- ggplot(model, aes(split_label, method_label, fill = displayed_rmse)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = fifelse(is.finite(displayed_rmse),
                                sprintf("%.3f", displayed_rmse), "Abstain")),
            size = 2.2) +
  scale_fill_gradient(low = "#F0F5F8", high = "#89ADC7", na.value = "#EEEEEE") +
  scale_x_discrete(labels = c("Cell line" = "Cell\nline", "Target state" = "Target\nstate",
    "Perturbation" = "Perturbation")) +
  labs(x = "Held-out unit", y = NULL, fill = "RMSE") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

p3b <- ggplot(model[
  method %in% c(
    "local_response_ridge", "kernel_transition",
    "conditional_rate_network"
  ) & is.finite(target_state_rmse)
], aes(coverage, target_state_rmse, colour = method, shape = split_label)) +
  geom_point(size = 2) +
  scale_colour_manual(values = method_colors, labels = c(
    local_response_ridge = "Local ridge", kernel_transition = "Kernel transition",
    conditional_rate_network = "Conditional rate")) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    title = "Abstention exposes the accuracy-coverage trade-off",
    x = "Non-abstained coverage",
    y = "Target-state RMSE",
    colour = NULL, shape = "Held-out unit"
  ) +
  theme(legend.position = "right")

p3c <- ggplot(model[coverage > 0 & is.finite(binned_probability_error)],
  aes(binned_probability_error, method_label, colour = split_label, shape = split_label)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Cell line" = "#4C78A8",
    "Target state" = "#E3A857", "Perturbation" = "#777777")) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15)) +
  labs(x = "Binned probability error", y = NULL, colour = NULL, shape = NULL) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

combo_metric <- fread("results/grouped_validation/norman_metrics.tsv")
combo_metric[, split_label := fifelse(split == "leave_gene", "Leave gene", "Leave family")]
p3d <- ggplot(combo_metric, aes(synergy_auc, method, colour = method, shape = split_label)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#A0A0A0") +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = method_colors) +
  scale_y_discrete(labels = c(nonlinear_xgboost = "XGBoost",
    kernel_regression = "Kernel regression", elastic_net = "Elastic net",
    additive_null = "Additive null")) +
  scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8)) +
  labs(x = "Synergy AUROC", y = NULL, shape = NULL) +
  guides(colour = "none") +
  theme(legend.position = "bottom")

fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(model, file.path(source_dir, "Fig3abc_model_metrics.tsv"), sep = "\t")
fwrite(combo_metric, file.path(source_dir, "Fig3d_combo_validation.tsv"),
       sep = "\t")
save_pub_r(fig3, "figures/main/Fig3_model_benchmark",
           width_mm = 183, height_mm = 145)

# Figure 4 -----------------------------------------------------------------
combo <- fread("results/norman_control/combination_interactions.tsv")
class_count <- combo[, .N, by = interaction_class]
class_count[, interaction_class := factor(
  interaction_class,
  levels = c("synergistic", "antagonistic", "redundant", "context-dependent")
)]
top_combo <- combo[order(-abs(replicate_projection_mean)),
                   head(.SD, 4L), by = interaction_class]
top_combo <- top_combo[order(replicate_projection_mean)]
top_combo[, combination := factor(
  combination, levels = rev(combination)
)]

p4a <- ggplot(class_count, aes(N, interaction_class,
                               fill = interaction_class)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = N), hjust = -0.3, size = 2.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "Most observed pairs deviate from simple additivity",
    x = "Combination conditions", y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

p4b <- ggplot(combo, aes(additive_movement, target_movement,
                         colour = interaction_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = palette_control[["neutral_mid"]]) +
  geom_point(alpha = 0.8, size = 1.25) +
  scale_colour_manual(values = class_colors) +
  coord_equal() +
  labs(
    title = "Observed target movement differs from the additive null",
    x = "Additive single-perturbation movement",
    y = "Observed pair movement",
    colour = NULL
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.2),
    legend.key.width = grid::unit(3, "mm")
  )

p4c <- ggplot(top_combo, aes(replicate_projection_mean, combination,
                             colour = interaction_class)) +
  geom_vline(xintercept = 0, colour = palette_control[["neutral_light"]]) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y", width = 0, linewidth = 0.4
  ) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = class_colors) +
  labs(
    title = "Eight gem groups resolve interaction direction",
    x = "Interaction projection (mean ± 95% CI)",
    y = NULL, colour = NULL
  ) +
  theme(legend.position = "none")

p4d <- ggplot(combo, aes(interaction_class, off_target_interaction,
                         colour = interaction_class)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.12, seed = 19),
             alpha = 0.6, size = 0.8) +
  scale_colour_manual(values = class_colors) +
  labs(
    title = "Interaction and off-target deviation are distinct",
    x = NULL,
    y = "Off-target interaction norm",
    colour = NULL
  ) +
  scale_x_discrete(labels = c(antagonistic = "Antagonistic", redundant = "Redundant",
    synergistic = "Synergistic", `context-dependent` = "Context-\ndependent")) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6))

fig4 <- p4a + p4b + p4c + p4d +
  plot_layout(design = "AB\nCB\nCD", heights = c(0.7, 0.7, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(class_count, file.path(source_dir, "Fig4a_interaction_counts.tsv"),
       sep = "\t")
fwrite(combo, file.path(source_dir, "Fig4bcd_combination_metrics.tsv"),
       sep = "\t")
fwrite(top_combo, file.path(source_dir, "Fig4c_displayed_pairs.tsv"), sep = "\t")
save_pub_r(fig4, "figures/main/Fig4_combination_control",
           width_mm = 183, height_mm = 150)

# Figure 5 -----------------------------------------------------------------
guide <- fread("results/multiperturb_multiome/target_control_roles.tsv")
guide <- guide[passes_task_qc == TRUE]
role_colors <- c(
  state_mover = palette_control[["genetic"]],
  barrier_breaker_candidate = palette_control[["dynamics"]],
  coupled_mover = palette_control[["multiome"]],
  weak_or_contextual = palette_control[["neutral_mid"]]
)
repro <- fread("results/multiperturb_multiome/target_reproducibility.tsv")
repro <- repro[passes_task_qc == TRUE]
p5a <- ggplot(guide, aes(rna_movement, atac_movement,
                         colour = control_role)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_colour_manual(values = role_colors) +
  labs(
    title = "Joint perturbation responses separate modality dominance",
    x = "RNA response movement",
    y = "Gene-linked ATAC response movement",
    colour = NULL
  ) +
  theme(legend.position = "right")

# One row per independently supported target; no duplicated movement dumbbell.
guide[, guide_gene := factor(guide_gene, levels = rev(sort(unique(guide_gene))))]
guide_cosine <- melt(guide[, .(guide_gene, median_rna_guide_cosine,
                               median_atac_guide_cosine)], id.vars = "guide_gene",
                     variable.name = "modality", value.name = "cosine")
p5b <- ggplot(guide_cosine, aes(modality, guide_gene, fill = cosine)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", cosine)), size = 2.2) +
  scale_x_discrete(labels = c("RNA", "ATAC")) +
  scale_fill_gradient(low = "#EDF4F4", high = "#77B3AE", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Guide cosine") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")
p5c <- ggplot(guide, aes(n_cells, guide_gene, colour = control_role)) +
  geom_segment(aes(x = 0, xend = n_cells, yend = guide_gene),
                colour = "#D9D9D9", linewidth = 0.5) +
  geom_point(aes(shape = factor(n_guides)), size = 2) +
  scale_colour_manual(values = role_colors,
    labels = c(coupled_mover = "Coupled", weak_or_contextual = "Weak / contextual")) +
  labs(x = "Cells per target", y = NULL, colour = NULL, shape = "Guides") +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0))
p5a <- p5a + theme(legend.position = "none")
fig5 <- p5a + p5b + p5c +
  plot_layout(widths = c(1.3, 1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(guide, file.path(source_dir, "Fig5ad_target_roles.tsv"), sep = "\t")
fwrite(repro, file.path(source_dir, "Fig5bc_guide_reproducibility.tsv"),
       sep = "\t")
save_pub_r(fig5, "figures/main/Fig5_multiome_control_roles",
           width_mm = 183, height_mm = 100)
fwrite(guide_cosine, file.path(source_dir, "Fig5b_guide_cosine_matrix.tsv"), sep = "\t")

# Figure 6 -----------------------------------------------------------------
deg <- fread("results/control_atlas/control_degeneracy.tsv")
deg[, landscape := fifelse(
  fragile_control_axis, "fragile",
  fifelse(redundant_control_landscape, "redundant", "intermediate")
)]
landscape_colors <- c(
  fragile = palette_control[["fail"]],
  intermediate = palette_control[["neutral_mid"]],
  redundant = palette_control[["multiome"]]
)
top_deg <- deg[order(-control_degeneracy)]
top_deg[, target_state_id := factor(
  target_state_id, levels = rev(target_state_id)
)]

p6a <- ggplot(
  deg[is.finite(pathway_substitutability)],
  aes(control_degeneracy, pathway_substitutability,
                       colour = landscape, size = observed_control_sets)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = landscape_colors) +
  scale_size_continuous(range = c(1, 4)) +
  labs(
    title = "Observed path diversity and disjoint-set similarity",
    x = "Effective number of response paths",
    y = "Disjoint-set endpoint similarity",
    colour = NULL, size = "Observed sets"
  ) +
  theme(legend.position = "right")

p6d <- ggplot(deg, aes(control_degeneracy, median_control_energy,
                       colour = landscape)) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = landscape_colors) +
  scale_y_log10() +
  labs(
    title = "Path diversity does not guarantee low relative effort",
    x = "Observed response-path diversity",
    y = "Median relative-effort proxy\n(log scale)",
    colour = NULL
  ) +
  theme(legend.position = "none")

sensitivity <- fread("results/control_atlas/control_degeneracy_sensitivity.tsv")
sensitivity[, basin := factor(sub("NORMAN_BASIN_", "", target_state_id),
  levels = sub("NORMAN_BASIN_", "", as.character(top_deg$target_state_id)))]
p6heat <- ggplot(sensitivity, aes(factor(cosine_threshold), basin,
                                 fill = control_degeneracy)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", control_degeneracy)), size = 2.1) +
  scale_fill_gradient(low = "#F1EDF6", high = "#9584B1") +
  labs(x = "Response-cosine threshold", y = "Response basin", fill = "Effective paths") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")
p6a <- p6a + scale_colour_manual(values = landscape_colors,
  labels = c(fragile = "Low diversity", intermediate = "Intermediate",
             redundant = "High diversity")) +
  guides(colour = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.box = "vertical")
fig6 <- p6heat + p6a + p6d +
  plot_layout(design = "AB\nAC", widths = c(1, 1.35)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(deg, file.path(source_dir, "Fig6_control_degeneracy.tsv"), sep = "\t")
fwrite(sensitivity, file.path(source_dir, "Fig6a_threshold_sensitivity.tsv"), sep = "\t")
save_pub_r(fig6, "figures/main/Fig6_control_degeneracy",
           width_mm = 183, height_mm = 145)

# Figure 7 -----------------------------------------------------------------
seqs <- fread("results/control_atlas/observed_support_sequences.tsv")
seqs[, early_max := pmax(p_apoptosis_6, p_proliferation_6, p_stress_6)]
seqs[, late_max := pmax(p_apoptosis_24, p_proliferation_24, p_stress_24)]
seqs[, supported := observed_support_policy ==
  "continue_same_intervention_to_24h"]
boundary <- fread("results/control_atlas/sequence_support_boundary.tsv")

paired <- melt(seqs[, .(cell_line, perturbation, supported, early_max, late_max)],
  id.vars = c("cell_line", "perturbation", "supported"),
  variable.name = "time", value.name = "probability")
paired[, time := factor(time, levels = c("early_max", "late_max"),
                        labels = c("6 h", "24 h"))]
p7a <- ggplot(paired, aes(time, probability,
    group = interaction(cell_line, perturbation), colour = supported)) +
  geom_line(alpha = 0.4, linewidth = 0.35) +
  geom_point(size = 0.5, alpha = 0.6) +
  facet_wrap(~perturbation, nrow = 1, strip.position = "bottom") +
  scale_colour_manual(values = c(`FALSE` = "#A0A0A0", `TRUE` = "#2A9D8F"),
                       labels = c("Abstain", "Supported continuation")) +
  labs(x = NULL, y = "Maximum program probability", colour = NULL) +
  theme(legend.position = "bottom")
p7b <- ggplot(seqs, aes(late_max - early_max, perturbation, colour = supported)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#A0A0A0") +
  geom_boxplot(aes(group = perturbation), colour = "#555555", outlier.shape = NA,
               width = 0.5, linewidth = 0.4) +
  geom_point(position = position_jitter(height = 0.12, width = 0, seed = 19),
             size = 1, alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "#A0A0A0", `TRUE` = "#2A9D8F")) +
  labs(x = "Change in maximum program probability (24 h - 6 h)", y = NULL) +
  theme(legend.position = "none")
boundary[, policy_label := c(same_drug_6_to_24h = "Same drug, 6 to 24 h",
  switch_drug_A_to_B = "Switch drug A to B", order_Norman_pair = "Order simultaneous pair",
  commitment_lock = "Commitment lock")[as.character(policy_space)]]
p7d <- ggplot(boundary, aes(status, policy_label)) +
  geom_point(aes(shape = status), colour = "#356C91", size = 3) +
  scale_shape_manual(values = c(abstain = 4, evaluated = 16)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")
fig7 <- p7a + p7b + p7d +
  plot_layout(design = "AA\nBC", heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(seqs, file.path(source_dir, "Fig7abc_sequences.tsv"), sep = "\t")
fwrite(paired, file.path(source_dir, "Fig7a_paired_timepoints.tsv"), sep = "\t")
fwrite(boundary, file.path(source_dir, "Fig7d_sequence_boundary.tsv"),
       sep = "\t")
save_pub_r(fig7, "figures/main/Fig7_observed_support_sequences",
           width_mm = 183, height_mm = 145)

# Figure 8 -----------------------------------------------------------------
task_summary <- fread("results/task_registry/task_gate_by_dataset.tsv")
task_summary[, dataset_label := c(
  mixseq_2020 = "MIX-Seq",
  norman_2019_combinatorial = "Norman",
  replogle_2022_k562_essential = "Replogle",
  yan_2024_multiperturbseq = "MultiPerturb"
)[dataset_id]]
evidence <- fread("results/control_atlas/control_node_evidence.tsv")
evidence_count <- evidence[, .N, by = evidence_class]
validation <- fread("results/grouped_validation/validation_coverage.tsv")
validation[, split := factor(split, levels = rev(split))]
validation[, status_label := fifelse(
  status == "complete", "complete", "type audit only"
)]
top_evidence <- evidence[
  causal_regulator == TRUE
][order(-datasets_with_causal_support, -synergistic_combinations,
        control_energy)][1:20]
top_evidence[, gene := factor(gene, levels = rev(gene))]

p8a <- ggplot(task_summary, aes(eligible_tasks,
                                reorder(dataset_label, eligible_tasks),
                                fill = dataset_label)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = eligible_tasks), hjust = -0.15, size = 2.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c(
    "Replogle" = palette_control[["genetic"]],
    "Norman" = palette_control[["combination"]],
    "MIX-Seq" = palette_control[["dynamics"]],
    "MultiPerturb" = palette_control[["multiome"]]
  )) +
  labs(
    title = sprintf(
      "The atlas contains %s quality-gated tasks",
      format(sum(task_summary$eligible_tasks), big.mark = ",")
    ),
    x = "Eligible state-transition tasks",
    y = NULL
  ) +
  theme(legend.position = "none")

p8b <- ggplot(evidence_count, aes(N, reorder(evidence_class, N),
                                  fill = evidence_class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = N), hjust = -0.15, size = 2.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(
    association_hub = palette_control[["neutral_mid"]],
    causal_regulator = palette_control[["genetic"]],
    state_marker = palette_control[["dynamics"]],
    validated_control_node_candidate = palette_control[["multiome"]],
    insufficient_causal_evidence = palette_control[["neutral_light"]]
  )) +
  labs(
    title = "Causal regulators are not automatically control nodes",
    x = "Candidate genes",
    y = NULL
  ) +
  theme(legend.position = "none")

# Explicit criteria distinguish breadth of causal support from control-node claims.
evidence_matrix <- melt(top_evidence[, .(gene,
  `Multiple datasets` = datasets_with_causal_support >= 2,
  `Multiple contexts` = independent_contexts_with_causal_support >= 2,
  `Path specific` = path_specific,
  `Pair consistent` = combination_consistent,
  `Control candidate` = validated_control_node_candidate)],
  id.vars = "gene", variable.name = "criterion", value.name = "supported")
p8c <- ggplot(evidence_matrix, aes(criterion, gene, fill = supported)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = c(`FALSE` = "#EEEEEE", `TRUE` = "#6C9FBF"),
                     labels = c("Not met", "Met")) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = NULL, fill = "Evidence criterion") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")
validation[, split_label := gsub("_", " ", sub("leave_", "", split))]
validation[, evidence_level := fifelse(status == "complete", "Numerical", "Audit only")]
p8d <- ggplot(validation, aes(evidence_level, split_label)) +
  geom_point(aes(shape = evidence_level), size = 2.6, colour = "#356C91") +
  scale_shape_manual(values = c(`Audit only` = 1, Numerical = 16)) +
  labs(x = "Validation evidence", y = "Held-out unit") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")
p8b <- p8b + scale_y_discrete(labels = c(causal_regulator = "Causal regulator",
  state_marker = "State marker", insufficient_causal_evidence = "Insufficient causal\nevidence"))
fig8 <- p8a + p8b + p8c + p8d +
  plot_layout(design = "AC\nBC\nDD", heights = c(1, 1, 1.1), widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
fwrite(task_summary, file.path(source_dir, "Fig8a_task_summary.tsv"),
       sep = "\t")
fwrite(evidence_count, file.path(source_dir, "Fig8b_evidence_classes.tsv"),
       sep = "\t")
fwrite(top_evidence, file.path(source_dir, "Fig8c_control_candidates.tsv"),
       sep = "\t")
fwrite(evidence_matrix, file.path(source_dir, "Fig8c_evidence_matrix.tsv"), sep = "\t")
fwrite(validation, file.path(source_dir, "Fig8d_validation_coverage.tsv"),
       sep = "\t")
save_pub_r(fig8, "figures/main/Fig8_atlas_boundaries",
           width_mm = 183, height_mm = 175)
