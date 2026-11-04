#!/usr/bin/env Rscript
# Plot combination responses, exposure contrasts and evidence availability.
# Run from the package root after run_all.R; ggplot2 and patchwork are optional.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
palette_control <- c(
  neutral_dark = "#2B2B2B",
  neutral_mid = "#7A7A7A",
  neutral_light = "#D9D9D9",
  genetic = "#4C78A8",
  combination = "#7B6BA8",
  multiome = "#2A9D8F",
  dynamics = "#E3A857",
  pass = "#2A9D8F",
  pending = "#A0A0A0",
  fail = "#C95A4B"
)

theme_nature_control <- function(base_size = 7, base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.4),
      legend.title = ggplot2::element_text(size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      strip.text = ggplot2::element_text(size = base_size, face = "plain"),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(7, 8, 7, 7),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}


theme_set(theme_nature_control())
source_dir <- "results/figure_exports"
dir.create(source_dir, recursive=TRUE, showWarnings=FALSE)
save_pub_r <- function(plot,filename,width_mm,height_mm) {
  stem <- file.path(source_dir,basename(filename))
  if (Sys.info()[["sysname"]]=="Darwin") quartz(type="pdf",file=paste0(stem,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Helvetica") else cairo_pdf(paste0(stem,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Helvetica")
  print(plot);dev.off()
  if (requireNamespace("ragg",quietly=TRUE)) {ragg::agg_png(paste0(stem,"_preview.png"),width=width_mm,height=height_mm,units="mm",res=160);print(plot);dev.off()}
}
class_colors <- c(synergistic="#2A9D8F",antagonistic="#C95A4B",redundant="#7B6BA8",`context-dependent`="#E3A857")
class_labels <- c(synergistic="Positive non-additive",antagonistic="Antagonistic",redundant="Redundant",`context-dependent`="Residual")

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
  scale_fill_manual(values = class_colors, labels = class_labels) +
  labs(
    title = "Most observed pairs deviate from simple additivity",
    x = "Combination conditions", y = NULL
  ) +
  theme(
    axis.text.y = element_text(colour="red"),
    legend.position = "none",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

p4a <- p4a + scale_y_discrete(labels=c(synergistic="Positive\nnon-additive",antagonistic="Antagonistic",redundant="Redundant",`context-dependent`="Residual"))
p4b <- ggplot(combo, aes(additive_movement, target_movement,
                         colour = interaction_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = palette_control[["neutral_mid"]]) +
  geom_point(alpha = 0.8, size = 1.25) +
  scale_colour_manual(values = class_colors, labels = class_labels) +
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
    legend.text = element_text(size = 5.2,colour="red"),
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
  scale_colour_manual(values = class_colors, labels = class_labels) +
  labs(
    title = "Eight gem groups resolve interaction direction",
    x = "Interaction projection (mean ± 95% technical interval)",
    y = NULL, colour = NULL
  ) +
  theme(legend.position = "none",axis.title.x=element_text(colour="red"))

p4d <- ggplot(combo, aes(interaction_class, off_target_interaction,
                         colour = interaction_class)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.12, seed = 19),
             alpha = 0.6, size = 0.8) +
  scale_colour_manual(values = class_colors, labels = class_labels) +
  labs(
    title = "Interaction and off-target deviation are distinct",
    x = NULL,
    y = "Off-target interaction norm",
    colour = NULL
  ) +
  scale_x_discrete(labels = c(antagonistic = "Antagonistic", redundant = "Redundant",
    synergistic = "Positive\nnon-additive", `context-dependent` = "Residual")) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6,colour="red"))

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
                       labels = c("Maximum did not increase", "Maximum increased")) +
  labs(x = NULL, y = "Maximum program probability", colour = NULL) +
  theme(legend.position = "bottom",legend.text=element_text(colour="red"))
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
  scale_x_discrete(labels=c(abstain="Not measured",evaluated="Duration contrast")) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x=element_text(colour="red"),axis.line = element_blank(), axis.ticks = element_blank(),
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
    x = "Eligible endpoint tasks",
    y = NULL
  ) +
  theme(legend.position = "none",axis.title.x=element_text(colour="red"))

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
    x = "Gene symbols",
    y = NULL
  ) +
  theme(legend.position = "none")

# Preserve the original displayed genes, while separating absent measurements
# and insufficient QC support from an evaluated criterion that was not met.
audit <- fread("results/multiome_validation/gene_criterion_four_state.tsv")
v <- merge(as.data.table(top_evidence),audit[,.(gene,multiple_datasets_state,multiple_contexts_state,specificity_state,combination_state,candidate_state)],by="gene",sort=FALSE)
v[,gene:=factor(gene,levels=levels(top_evidence$gene))]
evidence_matrix <- melt(v[,.(gene,`Multiple datasets`=multiple_datasets_state,`Multiple contexts`=multiple_contexts_state,`Path specific`=specificity_state,`Pair consistent`=combination_state,`Control candidate`=candidate_state)],id.vars="gene",variable.name="criterion",value.name="status")
evidence_matrix[status=="Not met",status:="Evaluated: not met"]
evidence_matrix[status=="Insufficient support",status:="Support insufficient"]
evidence_matrix[,status:=factor(status,levels=c("Met","Evaluated: not met","Support insufficient","Not measured"))]
p8c <- ggplot(evidence_matrix,aes(criterion,gene,fill=status))+
  geom_tile(colour="#BBBBBB",linewidth=.25,show.legend=TRUE)+
  scale_fill_manual(values=c("Met"="#6C9FBF","Evaluated: not met"="#D9A16E","Support insufficient"="#B8B8B8","Not measured"="#FFFFFF"),drop=FALSE)+
  scale_x_discrete(labels=function(x)gsub(" ","\n",x))+
  labs(x=NULL,y=NULL,fill=NULL)+guides(fill=guide_legend(nrow=2,byrow=TRUE,override.aes=list(fill=c("#6C9FBF","#D9A16E","#B8B8B8","#FFFFFF"))))+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom",legend.text=element_text(size=5.2,colour="red"),legend.key.size=grid::unit(2.5,"mm"))
validation[, split_label := gsub("_", " ", sub("leave_", "", split))]
validation[, evidence_level := fifelse(status == "complete", "Numerical", "Audit only")]
p8d <- ggplot(validation, aes(evidence_level, split_label)) +
  geom_point(aes(shape = evidence_level), size = 2.6, colour = "#356C91") +
  scale_shape_manual(values = c(`Audit only` = 1, Numerical = 16)) +
  labs(x = "Validation evidence", y = "Held-out unit") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")
p8b <- p8b + theme(axis.text.y=element_text(colour="red"),axis.title.x=element_text(colour="red")) + scale_y_discrete(labels = c(causal_regulator = "Operational\ncausal regulator",
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
