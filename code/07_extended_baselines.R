#!/usr/bin/env Rscript

# Extended baseline contract
# These baselines use the same Norman atomic splits and the same observed
# interaction endpoint as the fitted pair models. All predictors are constructed
# only from matched single-perturbation responses or public regulatory priors.
# Architecture-specific names are reserved for the original implementations.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})
set.seed(1908)
out_dir <- "results/extended_baselines"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_cosine <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(NA_real_)
  denom <- sqrt(sum(x[keep]^2) * sum(y[keep]^2))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(x[keep] * y[keep]) / denom
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = "spearman")
}

rank_auc <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- y[keep]
  score <- score[keep]
  if (!any(y) || all(y)) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y]) - sum(seq_len(sum(y)))) /
    (sum(y) * sum(!y))
}

best_threshold <- function(score, truth) {
  keep <- is.finite(score) & !is.na(truth)
  score <- score[keep]
  truth <- truth[keep]
  if (!length(score)) return(0)
  candidates <- unique(as.numeric(quantile(
    score, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE
  )))
  accuracy <- vapply(candidates, function(threshold) {
    mean((score > threshold) == truth)
  }, numeric(1))
  candidates[which.max(accuracy)]
}

combo <- fread("results/norman_control/combination_interactions.tsv")
condition_index <- fread("results/norman_control/condition_index.tsv")
response <- as.data.table(readRDS("results/norman_control/condition_response_delta.rds"))
response_genes <- setdiff(names(response), "condition_id")
response_matrix <- as.matrix(response[, ..response_genes])
rownames(response_matrix) <- response$condition_id

single_index <- condition_index[
  condition_type == "single" & !is.na(single_gene),
  .(condition_id, gene = single_gene)
]
single_index <- single_index[condition_id %in% rownames(response_matrix)]
single_gene_response <- rbindlist(lapply(
  split(single_index$condition_id, single_index$gene),
  function(ids) {
    values <- colMeans(response_matrix[ids, , drop = FALSE], na.rm = TRUE)
    as.data.table(as.list(values))
  }
), idcol = "gene")
single_genes <- single_gene_response$gene
single_matrix <- as.matrix(single_gene_response[, ..response_genes])
rownames(single_matrix) <- single_genes
single_matrix[!is.finite(single_matrix)] <- 0

pair_gene_set <- sort(unique(c(combo$gene_a, combo$gene_b)))
missing_single <- setdiff(pair_gene_set, rownames(single_matrix))
if (length(missing_single)) {
  stop("Missing matched single responses for: ",
       paste(missing_single, collapse = ", "))
}

single_norm <- sqrt(rowMeans(single_matrix^2))
single_scaled <- single_matrix / pmax(sqrt(rowSums(single_matrix^2)), 1e-12)
single_similarity <- tcrossprod(single_scaled)
diag(single_similarity) <- 1

# A response-similarity network is built without pair outcomes. The primary
# threshold is the 90th percentile of absolute off-diagonal similarity.
off_diagonal <- abs(single_similarity[upper.tri(single_similarity)])
network_threshold <- as.numeric(quantile(
  off_diagonal[is.finite(off_diagonal)], 0.90, na.rm = TRUE
))
adjacency <- abs(single_similarity) >= network_threshold
diag(adjacency) <- FALSE
network <- graph_from_adjacency_matrix(
  adjacency, mode = "undirected", diag = FALSE
)
edge_index <- which(upper.tri(adjacency) & adjacency, arr.ind = TRUE)
if (nrow(edge_index)) {
  E(network)$weight <- abs(single_similarity[edge_index])
  E(network)$distance <- 1 / pmax(E(network)$weight, 1e-6)
}
centrality <- strength(network, weights = E(network)$weight)
names(centrality) <- rownames(single_similarity)
distance_matrix <- distances(network, weights = E(network)$distance)
finite_distance <- distance_matrix[is.finite(distance_matrix) &
                                     distance_matrix > 0]
disconnected_distance <- if (length(finite_distance)) {
  max(finite_distance) + 1
} else {
  1
}

# Public CollecTRI interactions provide a regulator-activity ranking baseline.
# This is deliberately not called SCENIC because no fitted SCENIC regulon object
# exists for this dataset.
collectri_path <- "inputs/auxiliary.rds"
collectri_url <- "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&format=tsv"
collectri <- as.data.table(readRDS("inputs/auxiliary.rds")$collectri)
required_collectri <- c(
  "source_genesymbol", "target_genesymbol",
  "is_stimulation", "is_inhibition"
)
if (!all(required_collectri %in% names(collectri))) {
  stop("Unexpected CollecTRI schema.")
}
collectri <- unique(collectri[
  target_genesymbol %in% response_genes &
    !is.na(source_genesymbol) & nzchar(source_genesymbol),
  .(
    regulator = source_genesymbol,
    target = target_genesymbol,
    weight = fifelse(
      toupper(as.character(is_stimulation)) == "TRUE" &
        toupper(as.character(is_inhibition)) != "TRUE",
      1,
      fifelse(
        toupper(as.character(is_inhibition)) == "TRUE" &
          toupper(as.character(is_stimulation)) != "TRUE",
        -1, 0
      )
    )
  )
])
collectri <- collectri[, if (.N >= 5L) .SD, by = regulator]
regulators <- sort(unique(collectri$regulator))
regulon_activity <- matrix(
  0, nrow = nrow(single_matrix), ncol = length(regulators),
  dimnames = list(rownames(single_matrix), regulators)
)
for (reg_name in regulators) {
  edges <- collectri[regulator == reg_name]
  target_index <- match(edges$target, response_genes)
  regulon_activity[, reg_name] <- as.numeric(
    single_matrix[, target_index, drop = FALSE] %*% edges$weight
  ) / sqrt(nrow(edges))
}
regulon_norm <- sqrt(rowMeans(regulon_activity^2))
regulon_scaled <- regulon_activity /
  pmax(sqrt(rowSums(regulon_activity^2)), 1e-12)
regulon_similarity <- tcrossprod(regulon_scaled)

lookup_pair <- function(matrix, gene_a, gene_b) {
  matrix[cbind(match(gene_a, rownames(matrix)),
               match(gene_b, colnames(matrix)))]
}
combo[, `:=`(
  de_magnitude_score =
    single_norm[gene_a] + single_norm[gene_b],
  response_similarity_score =
    lookup_pair(single_similarity, gene_a, gene_b),
  response_network_centrality_score =
    centrality[gene_a] + centrality[gene_b],
  static_shortest_path_score = -vapply(
    seq_len(.N),
    function(i) {
      value <- distance_matrix[gene_a[i], gene_b[i]]
      if (is.finite(value)) value else disconnected_distance
    },
    numeric(1)
  ),
  collectri_regulon_score =
    regulon_norm[gene_a] + regulon_norm[gene_b] +
      lookup_pair(regulon_similarity, gene_a, gene_b)
)]

baseline_features <- c(
  differential_expression_ranking = "de_magnitude_score",
  response_similarity_ranking = "response_similarity_score",
  response_network_centrality = "response_network_centrality_score",
  static_shortest_path = "static_shortest_path_score",
  collectri_regulon_ranking = "collectri_regulon_score"
)

gene_counts <- table(c(combo$gene_a, combo$gene_b))
gene_order <- names(sort(gene_counts, decreasing = TRUE))
fold_load <- numeric(5)
gene_fold <- integer(length(gene_order))
names(gene_fold) <- gene_order
for (gene in gene_order) {
  fold <- which.min(fold_load)
  gene_fold[gene] <- fold
  fold_load[fold] <- fold_load[fold] + gene_counts[[gene]]
}
combo[, anchor_gene := fifelse(
  as.numeric(gene_counts[gene_a]) <= as.numeric(gene_counts[gene_b]),
  gene_a, gene_b
)]
combo[, evaluation_fold := unname(gene_fold[anchor_gene])]

family_map <- fread("results/grouped_validation/gene_family_map.tsv")
family_lookup <- setNames(family_map$primary_pfam, family_map$gene)
combo[, `:=`(
  family_a = unname(family_lookup[gene_a]),
  family_b = unname(family_lookup[gene_b])
)]
family_pair_count <- rbindlist(list(
  combo[family_a != "PFAM_UNANNOTATED",
        .(combination, family = family_a)],
  combo[family_b != "PFAM_UNANNOTATED",
        .(combination, family = family_b)]
))[, .(affected_pairs = uniqueN(combination)), by = family]
eligible_families <- family_pair_count[affected_pairs >= 3L, family]
family_count <- setNames(
  family_pair_count$affected_pairs, family_pair_count$family
)
combo[, family_anchor := vapply(seq_len(.N), function(i) {
  candidates <- unique(c(family_a[i], family_b[i]))
  candidates <- candidates[candidates %in% eligible_families]
  if (!length(candidates)) return(NA_character_)
  candidates[order(as.numeric(family_count[candidates]), candidates)][1]
}, character(1))]

fit_rank_baselines <- function(train, test, split, fold, held_out) {
  rbindlist(lapply(names(baseline_features), function(method) {
    feature <- baseline_features[[method]]
    train_score <- train[[feature]]
    test_score <- test[[feature]]
    keep <- is.finite(train_score) &
      is.finite(train$interaction_projection)
    fit <- lm(
      interaction_projection ~ score,
      data = data.frame(
        interaction_projection = train$interaction_projection[keep],
        score = train_score[keep]
      )
    )
    predicted <- as.numeric(predict(
      fit, newdata = data.frame(score = test_score)
    ))
    threshold <- best_threshold(
      as.numeric(predict(fit)),
      train$minimum_intervention_size[keep] == 2L
    )
    data.table(
      split, fold, held_out,
      combination = test$combination,
      method,
      feature,
      raw_score = test_score,
      predicted,
      decision_threshold = threshold,
      observed = test$interaction_projection,
      observed_class = test$interaction_class,
      observed_binary = test$interaction_class == "synergistic",
      observed_minimum_size = test$minimum_intervention_size,
      predicted_minimum_size = fifelse(predicted > threshold, 2L, 1L)
    )
  }))
}

predictions <- list()
for (fold in 1:5) {
  held_genes <- names(gene_fold)[gene_fold == fold]
  train <- combo[!gene_a %in% held_genes & !gene_b %in% held_genes]
  test <- combo[evaluation_fold == fold]
  if (nrow(train) >= 30L && nrow(test)) {
    predictions[[length(predictions) + 1L]] <- fit_rank_baselines(
      train, test, "leave_gene", fold, test$anchor_gene
    )
  }
}
families <- sort(unique(na.omit(combo$family_anchor)))
for (family in families) {
  held_genes <- family_map[primary_pfam == family, gene]
  train <- combo[!gene_a %in% held_genes & !gene_b %in% held_genes]
  test <- combo[family_anchor == family]
  if (nrow(train) >= 30L && nrow(test) >= 3L) {
    predictions[[length(predictions) + 1L]] <- fit_rank_baselines(
      train, test, "leave_gene_family", match(family, families),
      rep(family, nrow(test))
    )
  }
}
prediction <- rbindlist(predictions, use.names = TRUE, fill = TRUE)
if (anyDuplicated(prediction[, .(split, method, combination)])) {
  stop("Extended baseline predictions are not atomic and unique.")
}
fwrite(
  prediction,
  file.path(out_dir, "norman_baseline_predictions.tsv"),
  sep = "\t"
)

metrics <- prediction[, .(
  n = .N,
  interaction_rmse = sqrt(mean((predicted - observed)^2)),
  interaction_spearman = safe_spearman(predicted, observed),
  synergy_auc = rank_auc(observed_binary, predicted),
  constructed_intervention_size_accuracy = mean(
    predicted_minimum_size == observed_minimum_size
  )
), by = .(split, method)]
fwrite(
  metrics,
  file.path(out_dir, "norman_baseline_metrics.tsv"),
  sep = "\t"
)

existing_metrics <- fread(
  "results/grouped_validation/norman_metrics.tsv"
)
all_metrics <- rbindlist(
  list(existing_metrics, metrics), use.names = TRUE, fill = TRUE
)
fwrite(
  all_metrics,
  file.path(out_dir, "norman_all_method_metrics.tsv"),
  sep = "\t"
)

resource_provenance <- data.table(
  resource = "CollecTRI",
  url = collectri_url,
  local_path = collectri_path,
  bytes = file.info(collectri_path)$size,
  md5 = unname(tools::md5sum(collectri_path)),
  use = "signed_regulon_activity_ranking_not_SCENIC"
)
fwrite(
  resource_provenance,
  file.path(out_dir, "resource_provenance.tsv"),
  sep = "\t"
)

baseline_model_spec <- data.table(
  method = names(baseline_features),
  feature = unname(baseline_features),
  feature_provenance = c(
    "matched_single_perturbation_response_magnitude",
    "matched_single_perturbation_response_cosine",
    "strength_in_train_outcome_independent_single_response_graph",
    "inverse_cosine_distance_in_train_outcome_independent_single_response_graph",
    "signed_CollecTRI_projection_of_matched_single_perturbation_responses"
  ),
  uses_observed_pair_response_for_feature = FALSE,
  calibration_uses_training_outcomes_only = TRUE,
  architecture_claim = c(
    "descriptive_DE_rank_baseline",
    "descriptive_response_similarity_baseline",
    "descriptive_network_centrality_baseline",
    "descriptive_static_shortest_path_baseline",
    "public_prior_regulon_rank_baseline_not_SCENIC"
  )
)
fwrite(
  baseline_model_spec,
  file.path(out_dir, "baseline_model_spec.tsv"),
  sep = "\t"
)

