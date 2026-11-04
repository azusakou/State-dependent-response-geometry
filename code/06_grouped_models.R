#!/usr/bin/env Rscript

# Validation contract
# Every split operates on complete perturbation conditions. A test perturbation
# gene is removed from every training combination in leave-gene analyses.

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
})
set.seed(1907)
out_dir <- "results/grouped_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combo <- fread("results/norman_control/combination_interactions.tsv")
genes <- sort(unique(c(combo$gene_a, combo$gene_b)))
gene_family_map<-as.data.table(readRDS("inputs/auxiliary.rds")$family_map)
fwrite(
  gene_family_map,
  file.path(out_dir, "gene_family_map.tsv"),
  sep = "\t"
)
family_lookup <- setNames(
  gene_family_map$primary_pfam, gene_family_map$gene
)
combo[, `:=`(
  family_a = unname(family_lookup[gene_a]),
  family_b = unname(family_lookup[gene_b])
)]
feature_cols <- c(
  "single_a_movement", "single_b_movement", "additive_movement"
)
combo[, single_balance := pmin(single_a_movement, single_b_movement) /
  pmax(pmax(single_a_movement, single_b_movement), 1e-8)]
feature_cols <- c(feature_cols, "single_balance")
fwrite(
  data.table(
    feature = feature_cols,
    provenance = "matched_single_perturbations_only",
    uses_observed_pair_response = FALSE
  ),
  file.path(out_dir, "norman_model_spec.tsv"),
  sep = "\t"
)

fit_methods <- function(train, test) {
  xtr <- as.matrix(train[, ..feature_cols])
  xte <- as.matrix(test[, ..feature_cols])
  mu <- colMeans(xtr)
  sig <- apply(xtr, 2, sd)
  sig[!is.finite(sig) | sig < 1e-8] <- 1
  xtr <- sweep(sweep(xtr, 2, mu, "-"), 2, sig, "/")
  xte <- sweep(sweep(xte, 2, mu, "-"), 2, sig, "/")
  y <- train$interaction_projection

  elastic <- glmnet::cv.glmnet(
    xtr, y, alpha = 0.5,
    nfolds = min(5L, max(3L, floor(nrow(train) / 10))),
    standardize = FALSE
  )
  pred_elastic <- as.numeric(predict(elastic, xte, s = "lambda.1se"))
  train_elastic <- as.numeric(predict(elastic, xtr, s = "lambda.1se"))

  dtrain <- xgb.DMatrix(xtr, label = y)
  fit_xgb <- xgb.train(
    params = list(
      objective = "reg:squarederror", max_depth = 2, eta = 0.05,
      subsample = 0.8, colsample_bytree = 0.8, nthread = 2
    ),
    data = dtrain, nrounds = 120, verbose = 0
  )
  pred_xgb <- predict(fit_xgb, xte)
  train_xgb <- predict(fit_xgb, xtr)

  pred_graph <- vapply(seq_len(nrow(xte)), function(i) {
    d2 <- rowSums((xtr - rep(xte[i, ], each = nrow(xtr)))^2)
    idx <- order(d2)[seq_len(min(30L, length(d2)))]
    h <- sqrt(max(d2[idx], 1e-8))
    w <- exp(-d2[idx] / (2 * h^2))
    weighted.mean(y[idx], w)
  }, numeric(1))

  best_threshold <- function(score, truth) {
    candidates <- unique(as.numeric(quantile(
      score, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE
    )))
    acc <- vapply(candidates, function(th) {
      mean((score > th) == truth)
    }, numeric(1))
    candidates[which.max(acc)]
  }
  truth_two <- train$minimum_intervention_size == 2L
  th_elastic <- best_threshold(train_elastic, truth_two)
  th_xgb <- best_threshold(train_xgb, truth_two)
  th_graph <- best_threshold(y, truth_two)

  rbindlist(list(
    data.table(method = "additive_null",
               predicted = rep(0, nrow(test)),
               decision_threshold = .Machine$double.xmax),
    data.table(method = "elastic_net", predicted = pred_elastic,
               decision_threshold = th_elastic),
    data.table(method = "nonlinear_xgboost", predicted = pred_xgb,
               decision_threshold = th_xgb),
    data.table(method = "kernel_regression", predicted = pred_graph,
               decision_threshold = th_graph)
  ), idcol = FALSE)[, row_local := rep(seq_len(nrow(test)), 4)]
}

gene_counts <- table(c(combo$gene_a, combo$gene_b))
gene_order <- names(sort(gene_counts, decreasing = TRUE))
fold_load <- numeric(5)
gene_fold <- integer(length(gene_order))
names(gene_fold) <- gene_order
for (g in gene_order) {
  f <- which.min(fold_load)
  gene_fold[g] <- f
  fold_load[f] <- fold_load[f] + gene_counts[[g]]
}
combo[, anchor_gene := fifelse(
  as.numeric(gene_counts[gene_a]) <= as.numeric(gene_counts[gene_b]),
  gene_a, gene_b
)]
combo[, evaluation_fold := unname(gene_fold[anchor_gene])]

predictions <- list()
split_membership <- list()
for (fold in 1:5) {
  held <- names(gene_fold)[gene_fold == fold]
  train <- combo[!gene_a %in% held & !gene_b %in% held]
  test <- combo[evaluation_fold == fold]
  if (nrow(train) < 30L || !nrow(test)) next
  split_membership[[length(split_membership)+1L]] <- rbindlist(list(data.table(split="leave_gene",fold=fold,combination=train$combination,membership="training"),data.table(split="leave_gene",fold=fold,combination=test$combination,membership="test")))
  pred <- fit_methods(train, test)
  pred <- merge(
    pred, test[, .(
      row_local = .I, combination, gene_a, gene_b,
      anchor_gene,
      observed = interaction_projection,
      observed_class = interaction_class,
      observed_minimum_size = minimum_intervention_size
    )],
    by = "row_local"
  )
  pred[, `:=`(
    split = "leave_gene",
    fold = fold,
    held_out = anchor_gene
  )]
  predictions[[length(predictions) + 1L]] <- pred
}

family_pair_count <- rbindlist(list(
  combo[family_a != "PFAM_UNANNOTATED", .(combination, family = family_a)],
  combo[family_b != "PFAM_UNANNOTATED", .(combination, family = family_b)]
), use.names = TRUE)[, .(affected_pairs = uniqueN(combination)), by = family]
eligible_families <- family_pair_count[affected_pairs >= 3L, family]
family_count_lookup <- setNames(
  family_pair_count$affected_pairs, family_pair_count$family
)
combo[, family_anchor := vapply(seq_len(.N), function(i) {
  candidates <- unique(c(family_a[i], family_b[i]))
  candidates <- candidates[candidates %in% eligible_families]
  if (!length(candidates)) return(NA_character_)
  candidates[order(
    as.numeric(family_count_lookup[candidates]), candidates
  )][1]
}, character(1))]
families <- sort(unique(na.omit(combo$family_anchor)))
for (fam in families) {
  held_genes <- gene_family_map[primary_pfam == fam, gene]
  train <- combo[!gene_a %in% held_genes & !gene_b %in% held_genes]
  test <- combo[family_anchor == fam]
  if (nrow(train) < 30L || nrow(test) < 3L) next
  split_membership[[length(split_membership)+1L]] <- rbindlist(list(data.table(split="leave_gene_family",fold=match(fam,families),combination=train$combination,membership="training"),data.table(split="leave_gene_family",fold=match(fam,families),combination=test$combination,membership="test")))
  pred <- fit_methods(train, test)
  pred <- merge(
    pred, test[, .(
      row_local = .I, combination, gene_a, gene_b,
      family_anchor,
      observed = interaction_projection,
      observed_class = interaction_class,
      observed_minimum_size = minimum_intervention_size
    )],
    by = "row_local"
  )
  pred[, `:=`(
    split = "leave_gene_family",
    fold = match(fam, families),
    held_out = family_anchor
  )]
  predictions[[length(predictions) + 1L]] <- pred
}

fwrite(rbindlist(split_membership),file.path(out_dir,"norman_pair_membership.tsv"),sep="\t")
combo_pred <- rbindlist(predictions, use.names = TRUE, fill = TRUE)
combo_pred[, predicted_class := fifelse(
  predicted > 0, "synergistic", "non_synergistic"
)]
combo_pred[, observed_binary := observed_class == "synergistic"]
combo_pred[, predicted_minimum_size := fifelse(
  predicted > decision_threshold, 2L, 1L
)]
fwrite(combo_pred, file.path(out_dir, "norman_predictions.tsv"), sep = "\t")

rank_auc <- function(y, score) {
  pos <- score[y]
  neg <- score[!y]
  if (!length(pos) || !length(neg)) return(NA_real_)
  ranks <- rank(c(pos, neg))
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}

combo_metrics <- combo_pred[, .(
  n = .N,
  interaction_rmse = sqrt(mean((predicted - observed)^2)),
  interaction_spearman = safe_spearman(predicted, observed),
  synergy_auc = rank_auc(observed_binary, predicted),
  constructed_intervention_size_accuracy = mean(
    predicted_minimum_size == observed_minimum_size
  )
), by = .(split, method)]
fwrite(combo_metrics, file.path(out_dir, "norman_metrics.tsv"), sep = "\t")

# Leave-lineage validation on unique clone/future-condition units.
landscape <- as.data.table(readRDS("inputs/auxiliary.rds")$legacy_lineage)
state_features <- c(
  "cell_cycle_state", "stress_state", "chromatin_regulator_proxy",
  "metabolic_state", "mitochondrial_fraction"
)
lineage_predictions <- list()
for (dataset_nm in unique(landscape$dataset_id)) {
  dt <- landscape[dataset_id == dataset_nm]
  if (anyDuplicated(dt$analysis_unit_id)) {
    stop("Duplicate LARRY analysis units in ", dataset_nm)
  }
  lineages <- dt[, .N, by = dominant_lineage][N >= 30, dominant_lineage]
  for (lin in lineages) {
    test <- dt[dominant_lineage == lin]
    held_clones <- unique(test$clone_id)
    train <- dt[
      dominant_lineage != lin & !clone_id %in% held_clones
    ]
    xtr <- as.matrix(train[, ..state_features])
    xte <- as.matrix(test[, ..state_features])
    mu <- colMeans(xtr, na.rm = TRUE)
    sig <- apply(xtr, 2, sd, na.rm = TRUE)
    sig[!is.finite(sig) | sig < 1e-8] <- 1
    xtr <- sweep(sweep(xtr, 2, mu, "-"), 2, sig, "/")
    xte <- sweep(sweep(xte, 2, mu, "-"), 2, sig, "/")
    xtr[!is.finite(xtr)] <- 0
    xte[!is.finite(xte)] <- 0
    beta <- solve(
      crossprod(cbind(1, xtr)) +
        diag(c(0, rep(10, ncol(xtr)))),
      crossprod(cbind(1, xtr), train$p_differentiation)
    )
    pred <- as.numeric(cbind(1, xte) %*% beta)
    lineage_predictions[[length(lineage_predictions) + 1L]] <- data.table(
      dataset_id = dataset_nm,
      held_out_lineage = lin,
      clone_id = test$clone_id,
      analysis_unit_id = test$analysis_unit_id,
      observed = test$p_differentiation,
      predicted = pred
    )
  }
}
lineage_predictions <- rbindlist(lineage_predictions)
fwrite(lineage_predictions, file.path(out_dir, "leave_lineage_predictions.tsv"),
       sep = "\t")
lineage_metrics <- lineage_predictions[, .(
  n = .N,
  spearman = cor(observed, predicted, method = "spearman",
                 use = "complete.obs"),
  mse = mean((observed - predicted)^2, na.rm = TRUE)
), by = dataset_id]
fwrite(lineage_metrics, file.path(out_dir, "leave_lineage_metrics.tsv"),
       sep = "\t")

mix_metrics <- fread("results/mixseq_control_models/metrics.tsv")
larry_transfer <- as.data.table(readRDS("inputs/auxiliary.rds")$larry_transfer)

validation_coverage <- data.table(
  split = c(
    "leave_target_state", "leave_source_state", "leave_gene",
    "leave_gene_family", "leave_cell_line", "leave_lineage",
    "leave_perturbation_type", "leave_dataset"
  ),
  implementation = c(
    "MIXSeq grouped endpoint prediction",
    "MIXSeq grouped endpoint prediction",
    "Norman pair model; all held-out-gene pairs removed from training",
    paste(
      "Norman pair model; all genes sharing the held-out",
      "org.Hs.eg.db PFAM primary domain removed"
    ),
    "MIXSeq grouped endpoint prediction",
    "LARRY dominant descendant lineage",
    "atomic isolation audit only",
    "LARRY cytokine versus in-vitro transfer"
  ),
  status = c(
    "complete", "complete", "complete", "complete", "complete", "complete",
    "not_numerically_comparable", "complete"
  ),
  boundary = c(
    rep("complete perturbation condition is the atomic unit", 6),
    "genetic and chemical tasks lack a common target-state geometry",
    "complete clone/dataset units"
  )
)
fwrite(validation_coverage, file.path(out_dir, "validation_coverage.tsv"),
       sep = "\t")

print(combo_metrics)
print(lineage_metrics)
print(validation_coverage)
