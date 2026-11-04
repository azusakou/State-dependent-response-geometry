#!/usr/bin/env Rscript

# Model contract
# Local response ridge: neighbourhood-weighted ridge approximation.
# Kernel transition: kernel smoother over source-state neighbours.
# Conditional rate network: endpoint-derived conditional rate approximation.
# Unsupported perturbations and distant source states trigger abstention.

suppressPackageStartupMessages({
  library(data.table)
  library(nnet)
})
source("code/common.R")

set.seed(1904)
out_dir <- "results/mixseq_control_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- as.data.table(readRDS("inputs/mixseq_modeling.rds"))
task_registry <- fread("results/task_registry/state_transition_tasks.tsv")
eligible_ids <- task_registry[
  dataset_id == "mixseq_2020" & status == "eligible", task_id
]
dat <- dat[pair_id %in% eligible_ids]
dat[, primary_target := c("apoptosis", "proliferation", "stress")[
  max.col(as.matrix(.SD), ties.method = "first")
], .SDcols = c("p_apoptosis", "p_proliferation", "p_stress")]

z_cols <- grep("^z[0-9]+_mean$", names(dat), value = TRUE)
y_cols <- c("p_apoptosis", "p_proliferation", "p_stress")
pert_levels <- sort(unique(dat$perturbation))
# Source-state folds are learned from each matched-control latent profile, not
# from the cell-line label. This keeps leave-source-state distinct from
# leave-cell-line while retaining every complete perturbation condition.
baseline_z <- scale(as.matrix(dat[, ..z_cols]))
baseline_z[!is.finite(baseline_z)] <- 0
source_k <- min(5L, nrow(dat))
source_km <- kmeans(baseline_z, centers = source_k, nstart = 50)
dat[, source_state_cluster := source_km$cluster]

balanced_group_folds <- function(group, k = 5L) {
  tab <- sort(table(group), decreasing = TRUE)
  fold_load <- numeric(k)
  map <- integer(length(tab))
  names(map) <- names(tab)
  for (g in names(tab)) {
    f <- which.min(fold_load)
    map[g] <- f
    fold_load[f] <- fold_load[f] + tab[[g]]
  }
  unname(map[as.character(group)])
}

split_registry <- rbindlist(list(
  data.table(
    row_id = seq_len(nrow(dat)),
    split = "leave_cell_line",
    fold = balanced_group_folds(dat$cell_line, 5L)
  ),
  data.table(
    row_id = seq_len(nrow(dat)),
    split = "leave_source_state",
    fold = dat$source_state_cluster
  ),
  data.table(
    row_id = seq_len(nrow(dat)),
    split = "leave_perturbation",
    fold = match(dat$perturbation, pert_levels)
  ),
  data.table(
    row_id = seq_len(nrow(dat)),
    split = "leave_target_state",
    fold = match(dat$primary_target,
                 sort(unique(dat$primary_target)))
  )
))
fwrite(split_registry, file.path(out_dir, "split_assignments.tsv"), sep = "\t")

scale_train_test <- function(train, test) {
  mu <- colMeans(train)
  sig <- apply(train, 2, sd)
  sig[!is.finite(sig) | sig < 1e-8] <- 1
  list(
    train = sweep(sweep(train, 2, mu, "-"), 2, sig, "/"),
    test = sweep(sweep(test, 2, mu, "-"), 2, sig, "/")
  )
}

design_matrix <- function(d, levels) {
  z <- as.matrix(d[, ..z_cols])
  time <- log1p(fifelse(is.finite(d$time_hours), d$time_hours, 0))
  pert <- model.matrix(~ factor(perturbation, levels = levels) - 1, data = d)
  cbind(z, time = time, pert)
}

ridge_predict <- function(x_train, y_train, x_test, lambda = 1) {
  xx <- cbind(1, x_train)
  xt <- cbind(1, x_test)
  penalty <- diag(ncol(xx))
  penalty[1, 1] <- 0
  beta <- solve(crossprod(xx) + lambda * penalty,
                crossprod(xx, y_train))
  xt %*% beta
}

local_linear_predict <- function(z_train, x_train, y_train, z_test,
                                 x_test, k = 45L, lambda = 2) {
  pred <- matrix(NA_real_, nrow(z_test), ncol(y_train))
  unc <- numeric(nrow(z_test))
  for (i in seq_len(nrow(z_test))) {
    d2 <- rowSums((z_train - rep(z_test[i, ], each = nrow(z_train)))^2)
    idx <- order(d2)[seq_len(min(k, length(d2)))]
    h <- sqrt(max(d2[idx], 1e-8))
    w <- exp(-d2[idx] / (2 * h^2))
    xx <- cbind(1, x_train[idx, , drop = FALSE])
    penalty <- diag(ncol(xx))
    penalty[1, 1] <- 0
    beta <- solve(
      crossprod(xx, xx * w) + lambda * penalty,
      crossprod(xx, y_train[idx, , drop = FALSE] * w)
    )
    pred[i, ] <- c(1, x_test[i, ]) %*% beta
    residual <- y_train[idx, , drop = FALSE] - xx %*% beta
    unc[i] <- sqrt(mean(residual^2, na.rm = TRUE))
  }
  list(pred = pred, uncertainty = unc)
}

graph_predict <- function(z_train, y_train, z_test, train_pert, test_pert,
                          k = 60L) {
  pred <- matrix(NA_real_, nrow(z_test), ncol(y_train))
  unc <- numeric(nrow(z_test))
  support_n <- integer(nrow(z_test))
  for (i in seq_len(nrow(z_test))) {
    support <- which(train_pert == test_pert[i])
    support_n[i] <- length(support)
    if (!length(support)) {
      support <- seq_len(nrow(z_train))
    }
    d2 <- rowSums((z_train[support, , drop = FALSE] -
      rep(z_test[i, ], each = length(support)))^2)
    idx_local <- order(d2)[seq_len(min(k, length(d2)))]
    idx <- support[idx_local]
    h <- sqrt(max(d2[idx_local], 1e-8))
    w <- exp(-d2[idx_local] / (2 * h^2))
    w <- w / sum(w)
    pred[i, ] <- colSums(y_train[idx, , drop = FALSE] * w)
    unc[i] <- sqrt(sum(w * rowSums(
      (y_train[idx, , drop = FALSE] -
        rep(pred[i, ], each = length(idx)))^2
    )) / ncol(y_train))
  }
  list(pred = pred, uncertainty = unc, support_n = support_n)
}

neural_ode_predict <- function(x_train, y_train, time_train,
                               x_test, time_test, seeds = 1:3) {
  base_logit <- qlogis(0.1)
  tau_train <- pmax(time_train, 1) / 24
  tau_test <- pmax(time_test, 1) / 24
  rate <- sweep(qlogis(pmin(pmax(y_train, 1e-4), 1 - 1e-4)),
                2, base_logit, "-") /
    tau_train
  ensemble <- array(NA_real_, dim = c(nrow(x_test), ncol(y_train),
                                      length(seeds)))
  for (s in seq_along(seeds)) {
    set.seed(1904 + seeds[s])
    fit <- nnet::nnet(
      x_train, rate, size = 10, linout = TRUE, decay = 1e-2,
      maxit = 600, trace = FALSE, MaxNWts = 5000
    )
    rate_hat <- predict(fit, x_test)
    ensemble[, , s] <- plogis(base_logit + rate_hat * tau_test)
  }
  list(
    pred = apply(ensemble, c(1, 2), mean, na.rm = TRUE),
    uncertainty = apply(ensemble, 1, function(a) mean(apply(a, 1, sd),
                                                    na.rm = TRUE))
  )
}

nearest_support_threshold <- function(z) {
  if (nrow(z) < 3L) return(Inf)
  d <- as.matrix(dist(z))
  diag(d) <- Inf
  quantile(apply(d, 1, min), 0.95, na.rm = TRUE)
}

all_predictions <- list()
counter <- 0L
for (split_name in unique(split_registry$split)) {
  split_dt <- split_registry[split == split_name]
  for (fold_id in sort(unique(split_dt$fold))) {
    test_idx <- split_dt[fold == fold_id, row_id]
    train_idx <- setdiff(seq_len(nrow(dat)), test_idx)
    if (length(test_idx) == 0L || length(train_idx) < 20L) next

    x_raw_train <- design_matrix(dat[train_idx], pert_levels)
    x_raw_test <- design_matrix(dat[test_idx], pert_levels)
    x_scaled <- scale_train_test(x_raw_train, x_raw_test)
    z_scaled <- scale_train_test(
      as.matrix(dat[train_idx, ..z_cols]),
      as.matrix(dat[test_idx, ..z_cols])
    )
    y_train <- as.matrix(dat[train_idx, ..y_cols])
    y_test <- as.matrix(dat[test_idx, ..y_cols])
    time_train <- fifelse(is.finite(dat$time_hours[train_idx]),
                          dat$time_hours[train_idx], 24)
    time_test <- fifelse(is.finite(dat$time_hours[test_idx]),
                         dat$time_hours[test_idx], 24)

    unseen_pert <- !dat$perturbation[test_idx] %in%
      unique(dat$perturbation[train_idx])
    support_cut <- nearest_support_threshold(z_scaled$train)
    nearest_test <- vapply(seq_len(nrow(z_scaled$test)), function(i) {
      min(sqrt(rowSums((z_scaled$train -
        rep(z_scaled$test[i, ], each = nrow(z_scaled$train)))^2)))
    }, numeric(1))
    unsupported_state <- nearest_test > support_cut

    method_results <- list()
    global_mean <- matrix(
      rep(colMeans(y_train), each = length(test_idx)),
      nrow = length(test_idx)
    )
    method_results$training_mean <- list(
      pred = global_mean,
      uncertainty = rep(mean(apply(y_train, 2, sd)), length(test_idx))
    )

    ridge <- ridge_predict(
      x_scaled$train, y_train, x_scaled$test, lambda = 2
    )
    method_results$ridge_static <- list(
      pred = ridge,
      uncertainty = rep(sqrt(mean(
        (y_train - ridge_predict(
          x_scaled$train, y_train, x_scaled$train, lambda = 2
        ))^2
      )), length(test_idx))
    )

    elastic_pred <- matrix(NA_real_, nrow = length(test_idx), ncol = 3)
    for (j in seq_len(3)) {
      fit <- glmnet::cv.glmnet(
        x_scaled$train, y_train[, j], alpha = 0.5,
        nfolds = min(5L, max(3L, floor(length(train_idx) / 20))),
        standardize = FALSE
      )
      elastic_pred[, j] <- as.numeric(predict(
        fit, x_scaled$test, s = "lambda.1se"
      ))
    }
    method_results$elastic_net <- list(
      pred = elastic_pred,
      uncertainty = rep(NA_real_, length(test_idx))
    )

    method_results$local_response_ridge <- local_linear_predict(
      z_scaled$train, x_scaled$train, y_train,
      z_scaled$test, x_scaled$test
    )
    method_results$kernel_transition <- graph_predict(
      z_scaled$train, y_train, z_scaled$test,
      dat$perturbation[train_idx], dat$perturbation[test_idx]
    )
    method_results$conditional_rate_network <- neural_ode_predict(
      x_scaled$train, y_train, time_train,
      x_scaled$test, time_test
    )

    for (method in names(method_results)) {
      counter <- counter + 1L
      pred <- pmin(pmax(method_results[[method]]$pred, 0), 1)
      unc <- method_results[[method]]$uncertainty
      low_graph_support <- if (method == "kernel_transition") {
        method_results[[method]]$support_n < 10L
      } else {
        rep(FALSE, length(test_idx))
      }
      abstain <- if (method %in% c(
        "local_response_ridge", "kernel_transition",
        "conditional_rate_network"
      )) {
        unseen_pert | unsupported_state | low_graph_support |
          !is.finite(unc) | unc > 0.15
      } else {
        rep(FALSE, length(test_idx))
      }
      tmp <- data.table(
        split = split_name,
        fold = fold_id,
        method = method,
        row_id = test_idx,
        pair_id = dat$pair_id[test_idx],
        cell_line = dat$cell_line[test_idx],
        source_state_cluster = dat$source_state_cluster[test_idx],
        perturbation = dat$perturbation[test_idx],
        time_hours = dat$time_hours[test_idx],
        primary_target = dat$primary_target[test_idx],
        abstain = abstain,
        uncertainty = unc
      )
      for (j in seq_along(y_cols)) {
        nm <- sub("^p_", "", y_cols[j])
        tmp[, paste0("observed_", nm) := y_test[, j]]
        tmp[, paste0("predicted_", nm) := pred[, j]]
      }
      all_predictions[[counter]] <- tmp
    }
  }
}

predictions <- rbindlist(all_predictions, use.names = TRUE, fill = TRUE)
fwrite(predictions, file.path(out_dir, "predictions.tsv"), sep = "\t")

ece <- function(obs, pred, bins = 10L) {
  br <- seq(0, 1, length.out = bins + 1L)
  b <- cut(pred, br, include.lowest = TRUE, labels = FALSE)
  sum(vapply(seq_len(bins), function(k) {
    idx <- which(b == k)
    if (!length(idx)) return(0)
    length(idx) / length(obs) * abs(mean(obs[idx]) - mean(pred[idx]))
  }, numeric(1)))
}

metrics <- predictions[, {
  keep <- !abstain
  observed <- as.matrix(.SD[, paste0("observed_", c(
    "apoptosis", "proliferation", "stress"
  )), with = FALSE])
  predicted <- as.matrix(.SD[, paste0("predicted_", c(
    "apoptosis", "proliferation", "stress"
  )), with = FALSE])
  distance_obs <- sqrt(rowSums((observed - 0.1)^2))
  distance_pred <- sqrt(rowSums((predicted - 0.1)^2))
  primary_idx <- max.col(observed)
  off_obs <- rowSums(observed) - observed[cbind(seq_len(nrow(observed)),
                                                primary_idx)]
  off_pred <- rowSums(predicted) - predicted[cbind(seq_len(nrow(predicted)),
                                                   primary_idx)]
  rare_obs <- apply(observed, 1, max) >= 0.8
  rare_pred <- apply(predicted, 1, max) >= 0.8
  list(
    n = .N,
    evaluated_n = sum(keep),
    coverage = mean(keep),
    target_state_rmse = if (any(keep)) sqrt(mean(
      (predicted[keep, ] - observed[keep, ])^2
    )) else NA_real_,
    target_distance_error = if (any(keep)) mean(abs(
      distance_pred[keep] - distance_obs[keep]
    )) else NA_real_,
    relative_effort_proxy_error = if (any(keep)) mean(abs(
      1 / (distance_pred[keep]^2 + 1e-8) -
        1 / (distance_obs[keep]^2 + 1e-8)
    )) else NA_real_,
    non_target_program_divergence = if (any(keep)) mean(abs(
      off_pred[keep] - off_obs[keep]
    )) else NA_real_,
    target_program_threshold_accuracy = if (any(keep)) mean(
      (apply(predicted[keep, , drop = FALSE], 1, max) >= 0.5) ==
        (apply(observed[keep, , drop = FALSE], 1, max) >= 0.5)
    ) else NA_real_,
    high_program_probability_recall = if (sum(rare_obs & keep) > 0) {
      sum(rare_pred & rare_obs & keep) / sum(rare_obs & keep)
    } else NA_real_,
    binned_probability_error = if (any(keep)) mean(vapply(
      seq_len(3), function(j) ece(observed[keep, j], predicted[keep, j]),
      numeric(1)
    )) else NA_real_
  )
}, by = .(split, method)]
fwrite(metrics, file.path(out_dir, "metrics.tsv"), sep = "\t")
print(metrics)
