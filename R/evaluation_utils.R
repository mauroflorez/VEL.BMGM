#' Summary of VEL-BMGM Model Fit
#'
#' Prints a concise summary of the fitted model including predictor selection,
#' covariate selection, intercept, and graph structure.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold for selection (default 0.5).
#' @param credible_level Credible interval level (default 0.95).
#' @param predictor_names Optional character vector of predictor names.
#' @param covariate_names Optional character vector of covariate names.
#'
#' @return Invisibly returns a list with summary components.
#' @export
summary_bmgm <- function(fit, threshold = 0.5, credible_level = 0.95,
                         predictor_names = NULL, covariate_names = NULL) {

  nburn   <- fit$nburn
  nsample <- fit$nsample
  n       <- nrow(fit$X)
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  alpha   <- 1 - credible_level

  if (is.null(predictor_names))  predictor_names  <- paste0("x", 1:p)
  if (is.null(covariate_names))  covariate_names  <- paste0("z", 1:K)

  # --- Predictor PPI ---
  post_idx    <- (nburn + 1):(nburn + nsample)
  gamma_post  <- fit$post_gamma[post_idx, 1:p, drop = FALSE]
  ppi_pred    <- colMeans(gamma_post)

  # --- Covariate PPI ---
  gt_post <- fit$post_gamma_tilde[1:p, , post_idx, drop = FALSE]
  ppi_cov <- apply(gt_post, c(1, 2), mean)

  # --- Intercept ---
  beta0_post <- fit$post_beta0[post_idx]
  beta0_mean <- mean(beta0_post)
  beta0_ci   <- quantile(beta0_post, probs = c(alpha / 2, 1 - alpha / 2))

  # --- tau_w ---
  tauw_post <- fit$post_tauw[post_idx]
  tauw_mean <- mean(tauw_post)

  # --- Graph ---
  G_post   <- fit$post_G[post_idx, , drop = FALSE]
  edge_ppi <- colMeans(G_post)
  n_edges  <- sum(edge_ppi > threshold)
  q        <- ncol(fit$post_gamma)
  max_edges <- q * (q - 1) / 2

  # --- Print ---
  cat("=== VEL-BMGM Model Summary ===\n\n")
  cat(sprintf("  Observations: %d | Predictors: %d | Covariates: %d\n", n, p, K))
  cat(sprintf("  MCMC: %d burn-in + %d samples\n", nburn, nsample))
  cat(sprintf("  Intercept: %.3f  [%.3f, %.3f]\n", beta0_mean, beta0_ci[1], beta0_ci[2]))
  cat(sprintf("  tau_w (GP weight precision): %.3f\n\n", tauw_mean))

  cat("--- Predictor Selection (Level 1) ---\n")
  pred_df <- data.frame(
    Predictor = predictor_names,
    Type      = fit$type,
    PPI       = round(ppi_pred, 4),
    Selected  = ifelse(ppi_pred > threshold, "*", "")
  )
  print(pred_df, row.names = FALSE)

  selected_j <- which(ppi_pred > threshold)
  cat(sprintf("\n  %d / %d predictors selected (threshold = %.2f)\n\n", length(selected_j), p, threshold))

  if (length(selected_j) > 0) {
    cat("--- Covariate Selection (Level 2) ---\n")
    for (j in selected_j) {
      cat(sprintf("  %s:", predictor_names[j]))
      for (k in 1:K) {
        flag <- ifelse(ppi_cov[j, k] > threshold, "*", " ")
        cat(sprintf("  %s=%.3f%s", covariate_names[k], ppi_cov[j, k], flag))
      }
      cat("\n")
    }
    cat("\n")
  }

  cat(sprintf("--- Graph: %d / %d edges selected ---\n", n_edges, max_edges))

  out <- list(
    ppi_predictors = ppi_pred,
    ppi_covariates = ppi_cov,
    beta0_mean     = beta0_mean,
    beta0_ci       = beta0_ci,
    tauw_mean      = tauw_mean,
    n_edges        = n_edges,
    predictor_names = predictor_names,
    covariate_names = covariate_names,
    threshold       = threshold
  )
  invisible(out)
}


#' Posterior Inclusion Probability Table
#'
#' Produces a publication-ready data frame of PPIs at both selection levels.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold for selection (default 0.5).
#' @param predictor_names Optional character vector of predictor names.
#' @param covariate_names Optional character vector of covariate names.
#'
#' @return A data frame with predictor and covariate PPIs.
#' @export
pip_table <- function(fit, threshold = 0.5,
                      predictor_names = NULL, covariate_names = NULL) {

  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names))  predictor_names  <- paste0("x", 1:p)
  if (is.null(covariate_names))  covariate_names  <- paste0("z", 1:K)

  # Level-1 PPI
  gamma_post <- fit$post_gamma[post_idx, 1:p, drop = FALSE]
  ppi_pred   <- colMeans(gamma_post)

  # Level-2 PPI
  gt_post <- fit$post_gamma_tilde[1:p, , post_idx, drop = FALSE]
  ppi_cov <- apply(gt_post, c(1, 2), mean)

  df <- data.frame(
    Predictor      = predictor_names,
    Type           = fit$type,
    PPI_predictor  = round(ppi_pred, 4),
    Selected       = ppi_pred > threshold
  )

  for (k in 1:K) {
    df[[paste0("PPI_", covariate_names[k])]] <- round(ppi_cov[, k], 4)
  }

  return(df)
}


#' Plot Varying Coefficients
#'
#' Plots estimated varying-coefficient functions beta_j(Z) with credible bands.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param predictor_indices Which predictors to plot (default: all selected with PPI > 0.5).
#' @param covariate_index Which covariate dimension to use as x-axis.
#' @param credible_level Credible interval level (default 0.95).
#' @param predictor_names Optional character vector of predictor names.
#' @param covariate_names Optional character vector of covariate names.
#' @param true_functions Optional named list of functions for overlay (e.g., list("x1" = function(z) 2*z^2+0.5)).
#'
#' @return A ggplot object (one facet per predictor).
#' @export
plot_varying_coefficients <- function(fit, predictor_indices = NULL,
                                     covariate_index = 1,
                                     credible_level = 0.95,
                                     predictor_names = NULL,
                                     covariate_names = NULL,
                                     true_functions = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for this function")

  nburn   <- fit$nburn
  nsample <- fit$nsample
  n       <- nrow(fit$X)
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  alpha   <- 1 - credible_level
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names))  predictor_names  <- paste0("x", 1:p)
  if (is.null(covariate_names))  covariate_names  <- paste0("z", 1:K)

  # Auto-select predictors by PPI
  if (is.null(predictor_indices)) {
    ppi <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])
    predictor_indices <- which(ppi > 0.5)
    if (length(predictor_indices) == 0) {
      message("No predictors selected at PPI > 0.5. Using top 3 by PPI.")
      predictor_indices <- order(ppi, decreasing = TRUE)[1:min(3, p)]
    }
  }

  # Z values for sorting
  Z_orig <- fit$Z
  z_vals <- Z_orig[, covariate_index]
  sort_order <- order(z_vals)

  plot_data <- list()
  for (j in predictor_indices) {
    beta_samples <- fit$beta_j[, j, post_idx]  # n x nsample
    beta_mean <- rowMeans(beta_samples)
    beta_lo   <- apply(beta_samples, 1, quantile, probs = alpha / 2)
    beta_hi   <- apply(beta_samples, 1, quantile, probs = 1 - alpha / 2)

    df_j <- data.frame(
      z          = z_vals[sort_order],
      beta_mean  = beta_mean[sort_order],
      beta_lo    = beta_lo[sort_order],
      beta_hi    = beta_hi[sort_order],
      predictor  = predictor_names[j]
    )
    plot_data[[length(plot_data) + 1]] <- df_j
  }

  plot_df <- do.call(rbind, plot_data)

  g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = z)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = beta_lo, ymax = beta_hi), alpha = 0.3, fill = "steelblue") +
    ggplot2::geom_line(ggplot2::aes(y = beta_mean), color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::facet_wrap(~predictor, scales = "free_y") +
    ggplot2::labs(
      x = covariate_names[covariate_index],
      y = expression(hat(beta)[j](z)),
      title = "Estimated Varying Coefficients"
    ) +
    ggplot2::theme_minimal()

  # Overlay true functions if provided
  if (!is.null(true_functions)) {
    for (j in predictor_indices) {
      name_j <- predictor_names[j]
      if (name_j %in% names(true_functions)) {
        z_grid <- seq(min(z_vals), max(z_vals), length.out = 200)
        true_vals <- true_functions[[name_j]](z_grid)
        true_df <- data.frame(z = z_grid, y_true = true_vals, predictor = name_j)
        g <- g + ggplot2::geom_line(data = true_df,
                                    ggplot2::aes(x = z, y = y_true),
                                    color = "red", linetype = "dashed", linewidth = 0.7)
      }
    }
  }

  return(g)
}


#' Plot PPI Heatmap
#'
#' Heatmap of two-level variable selection results. Rows are predictors,
#' columns are covariates, color encodes covariate-level PPI.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold for annotation (default 0.5).
#' @param predictor_names Optional character vector of predictor names.
#' @param covariate_names Optional character vector of covariate names.
#'
#' @return A ggplot heatmap object.
#' @export
plot_pip_heatmap <- function(fit, threshold = 0.5,
                             predictor_names = NULL, covariate_names = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for this function")

  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names))  predictor_names  <- paste0("x", 1:p)
  if (is.null(covariate_names))  covariate_names  <- paste0("z", 1:K)

  # Level-1 PPI
  ppi_pred <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])

  # Level-2 PPI
  gt_post <- fit$post_gamma_tilde[1:p, , post_idx, drop = FALSE]
  ppi_cov <- apply(gt_post, c(1, 2), mean)

  # Build long-form data
  df <- expand.grid(
    Covariate = covariate_names,
    Predictor = predictor_names,
    stringsAsFactors = FALSE
  )
  df$PPI <- as.vector(t(ppi_cov))
  df$PPI_pred <- rep(ppi_pred, each = K)

  # Order predictors by Level-1 PPI (descending)
  pred_order <- predictor_names[order(ppi_pred, decreasing = TRUE)]
  df$Predictor <- factor(df$Predictor, levels = pred_order)
  df$Covariate <- factor(df$Covariate, levels = covariate_names)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = Covariate, y = Predictor, fill = PPI)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(PPI, 2)), size = 3) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "lightyellow",
                                  high = "steelblue", midpoint = 0.25,
                                  limits = c(0, 1)) +
    ggplot2::labs(
      title = "Covariate-level PPI Heatmap",
      subtitle = sprintf("Rows ordered by predictor PPI (threshold = %.2f)", threshold),
      fill = "PPI"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  return(g)
}


#' Plot Predictor Graph
#'
#' Visualize the estimated predictor graph as a network diagram.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold for edge inclusion (default 0.5).
#' @param predictor_names Optional character vector of predictor names.
#' @param node_color Optional color vector for nodes (length p).
#' @param layout igraph layout name (default "circle").
#'
#' @return An igraph plot (invisible igraph object).
#' @export
plot_graph <- function(fit, threshold = 0.5,
                       predictor_names = NULL,
                       node_color = NULL,
                       layout = "circle") {

  if (!requireNamespace("igraph", quietly = TRUE))
    stop("igraph is required for this function")

  adj <- fit$adj_G
  p   <- nrow(adj)

  if (is.null(predictor_names)) predictor_names <- paste0("x", 1:p)

  # Threshold the adjacency
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0

  colnames(adj_bin) <- predictor_names
  rownames(adj_bin) <- predictor_names

  g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)

  # Node color by PPI if not provided
  if (is.null(node_color)) {
    nburn   <- fit$nburn
    nsample <- fit$nsample
    post_idx <- (nburn + 1):(nburn + nsample)
    ppi <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])
    node_color <- ifelse(ppi > 0.5, "steelblue", "grey80")
  }

  layout_fn <- switch(layout,
    "circle" = igraph::layout_in_circle,
    "fr"     = igraph::layout_with_fr,
    "kk"     = igraph::layout_with_kk,
    igraph::layout_in_circle
  )

  igraph::plot.igraph(g,
    vertex.color = node_color,
    vertex.size  = 25,
    vertex.label = predictor_names,
    vertex.label.cex = 0.8,
    edge.width   = 2,
    layout       = layout_fn(g),
    main         = "Estimated Predictor Graph"
  )

  invisible(g)
}


#' Plot PPI Bar Chart
#'
#' Bar plot of Level-1 posterior inclusion probabilities with threshold line.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold (default 0.5).
#' @param predictor_names Optional character vector of predictor names.
#'
#' @return A ggplot object.
#' @export
plot_pip_bar <- function(fit, threshold = 0.5, predictor_names = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for this function")

  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(fit$X)
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names)) predictor_names <- paste0("x", 1:p)

  ppi <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])

  df <- data.frame(
    Predictor = factor(predictor_names, levels = predictor_names[order(ppi, decreasing = TRUE)]),
    PPI       = ppi,
    Selected  = ppi > threshold
  )

  g <- ggplot2::ggplot(df, ggplot2::aes(x = Predictor, y = PPI, fill = Selected)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70")) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "Predictor Posterior Inclusion Probabilities",
                  y = "PPI", x = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::guides(fill = "none")

  return(g)
}


#' MCMC Traceplots
#'
#' Produces traceplots for key model parameters.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param parameter One of "beta0", "gamma", "tauw", "r".
#' @param predictor_indices For gamma or r: which predictors (default: first 4).
#' @param include_burnin Whether to include burn-in iterations (default TRUE).
#'
#' @return A ggplot object.
#' @export
plot_trace <- function(fit, parameter = "beta0",
                       predictor_indices = NULL,
                       include_burnin = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for this function")

  nburn   <- fit$nburn
  nsample <- fit$nsample
  total   <- nburn + nsample
  p       <- ncol(fit$X)

  if (include_burnin) {
    iter_range <- 1:total
  } else {
    iter_range <- (nburn + 1):total
  }

  if (parameter == "beta0") {
    df <- data.frame(Iteration = iter_range, Value = fit$post_beta0[iter_range])
    g <- ggplot2::ggplot(df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(alpha = 0.5, linewidth = 0.3) +
      ggplot2::labs(title = expression("Traceplot: " * beta[0]), y = expression(beta[0]))

  } else if (parameter == "tauw") {
    df <- data.frame(Iteration = iter_range, Value = fit$post_tauw[iter_range])
    g <- ggplot2::ggplot(df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(alpha = 0.5, linewidth = 0.3) +
      ggplot2::labs(title = expression("Traceplot: " * tau[w]), y = expression(tau[w]))

  } else if (parameter == "gamma") {
    if (is.null(predictor_indices)) predictor_indices <- 1:min(4, p)
    dfs <- lapply(predictor_indices, function(j) {
      data.frame(Iteration = iter_range,
                 Value = fit$post_gamma[iter_range, j],
                 Predictor = paste0("x", j))
    })
    df <- do.call(rbind, dfs)
    g <- ggplot2::ggplot(df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(alpha = 0.4, linewidth = 0.3) +
      ggplot2::facet_wrap(~Predictor) +
      ggplot2::labs(title = expression("Traceplot: " * gamma[j]), y = expression(gamma[j]))

  } else if (parameter == "r") {
    if (is.null(predictor_indices)) predictor_indices <- 1:min(4, p)
    dfs <- lapply(predictor_indices, function(j) {
      data.frame(Iteration = iter_range,
                 Value = fit$post_r[iter_range, j],
                 Predictor = paste0("x", j))
    })
    df <- do.call(rbind, dfs)
    g <- ggplot2::ggplot(df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(alpha = 0.4, linewidth = 0.3) +
      ggplot2::facet_wrap(~Predictor, scales = "free_y") +
      ggplot2::labs(title = expression("Traceplot: " * r[j]), y = expression(r[j]))

  } else {
    stop("parameter must be one of: 'beta0', 'gamma', 'tauw', 'r'")
  }

  if (include_burnin) {
    g <- g + ggplot2::geom_vline(xintercept = nburn, linetype = "dashed",
                                  color = "red", alpha = 0.6)
  }

  g <- g + ggplot2::theme_minimal()
  return(g)
}


#' Classification Metrics
#'
#' Computes standard binary classification metrics from model predictions.
#'
#' @param y_true True binary labels.
#' @param y_pred_prob Predicted probabilities (from predict_jbmgm).
#' @param threshold Classification threshold (default 0.5).
#'
#' @return A named vector of metrics.
#' @export
classification_metrics <- function(y_true, y_pred_prob, threshold = 0.5) {
  y_pred <- as.integer(y_pred_prob > threshold)

  TP <- sum(y_pred == 1 & y_true == 1)
  TN <- sum(y_pred == 0 & y_true == 0)
  FP <- sum(y_pred == 1 & y_true == 0)
  FN <- sum(y_pred == 0 & y_true == 1)

  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  f1          <- ifelse((precision + sensitivity) > 0,
                        2 * precision * sensitivity / (precision + sensitivity), 0)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denom > 0, (TP * TN - FP * FN) / denom, 0)

  brier <- mean((y_pred_prob - y_true)^2)

  auc_val <- NA
  if (requireNamespace("pROC", quietly = TRUE)) {
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(y_true, y_pred_prob, quiet = TRUE))),
      error = function(e) NA
    )
  }

  c(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
    Precision = precision, F1 = f1, MCC = mcc, AUC = auc_val, Brier = brier)
}


#' Selection Metrics
#'
#' Computes predictor selection accuracy given true predictor indices.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param true_predictors Integer vector of true predictor indices.
#' @param threshold PPI threshold for selection (default 0.5).
#'
#' @return A named vector of selection metrics.
#' @export
selection_metrics <- function(fit, true_predictors, threshold = 0.5) {
  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(fit$X)
  post_idx <- (nburn + 1):(nburn + nsample)

  ppi      <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])
  selected <- as.integer(ppi > threshold)

  true_labels <- rep(0, p)
  true_labels[true_predictors] <- 1

  TP <- sum(selected == 1 & true_labels == 1)
  TN <- sum(selected == 0 & true_labels == 0)
  FP <- sum(selected == 1 & true_labels == 0)
  FN <- sum(selected == 0 & true_labels == 1)

  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  f1          <- ifelse((precision + sensitivity) > 0,
                        2 * precision * sensitivity / (precision + sensitivity), 0)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denom > 0, (TP * TN - FP * FN) / denom, 0)

  c(Sensitivity = sensitivity, Specificity = specificity,
    Precision = precision, F1 = f1, MCC = mcc)
}


#' Posterior Predictive Check
#'
#' Generates posterior predictive samples and computes Bayesian p-values.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param n_rep Number of posterior predictive replications (default 500).
#' @param summary_fn Summary statistic function (default: mean).
#'
#' @return A list with replicated stats, observed stat, and Bayesian p-value.
#' @export
posterior_predictive_check <- function(fit, n_rep = 500,
                                      summary_fn = mean) {
  nburn   <- fit$nburn
  nsample <- fit$nsample
  n       <- nrow(fit$X)
  p       <- ncol(fit$X)
  Y       <- fit$Y
  post_idx <- (nburn + 1):(nburn + nsample)

  sample_iters <- sample(post_idx, size = min(n_rep, nsample))

  T_rep <- numeric(length(sample_iters))
  for (s in seq_along(sample_iters)) {
    m <- sample_iters[s]
    active <- which(fit$post_gamma[m, 1:p] == 1)
    beta_0_s <- fit$post_beta0[m]

    eta_s <- rep(beta_0_s, n)
    if (length(active) > 0) {
      beta_m <- matrix(fit$beta_j[, active, m], nrow = n, ncol = length(active))
      eta_s <- eta_s + rowSums(fit$X[, active, drop = FALSE] * beta_m)
    }
    prob_s <- 1 / (1 + exp(-eta_s))
    y_rep  <- rbinom(n, 1, prob_s)
    T_rep[s] <- summary_fn(y_rep)
  }

  T_obs   <- summary_fn(Y)
  p_value <- mean(T_rep >= T_obs)

  list(T_rep = T_rep, T_obs = T_obs, bayesian_pvalue = p_value)
}


#' Coefficient Summary Table
#'
#' Summarizes varying coefficients beta_j(Z) — posterior means and credible intervals.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param predictor_indices Which predictors to summarize (default: selected).
#' @param credible_level Credible interval level (default 0.95).
#' @param predictor_names Optional character vector of predictor names.
#'
#' @return A data frame with pointwise summaries for each observation.
#' @export
coef_summary <- function(fit, predictor_indices = NULL,
                         credible_level = 0.95,
                         predictor_names = NULL) {

  nburn   <- fit$nburn
  nsample <- fit$nsample
  n       <- nrow(fit$X)
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  alpha   <- 1 - credible_level
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names)) predictor_names <- paste0("x", 1:p)

  if (is.null(predictor_indices)) {
    ppi <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])
    predictor_indices <- which(ppi > 0.5)
  }

  results <- list()
  for (j in predictor_indices) {
    beta_samples <- fit$beta_j[, j, post_idx]  # n x nsample
    df_j <- data.frame(
      obs        = 1:n,
      predictor  = predictor_names[j],
      post_mean  = rowMeans(beta_samples),
      post_sd    = apply(beta_samples, 1, sd),
      ci_lower   = apply(beta_samples, 1, quantile, probs = alpha / 2),
      ci_upper   = apply(beta_samples, 1, quantile, probs = 1 - alpha / 2)
    )
    # Add covariate values
    for (k in 1:K) {
      df_j[[paste0("Z", k)]] <- fit$Z[, k]
    }
    df_j$significant <- (df_j$ci_lower > 0) | (df_j$ci_upper < 0)
    results[[length(results) + 1]] <- df_j
  }

  do.call(rbind, results)
}


#' Edge Posterior Inclusion Probability Table
#'
#' Produces a table of edge PPIs for the predictor graph.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param threshold PPI threshold for edge inclusion (default 0.5).
#' @param predictor_names Optional character vector of predictor names.
#'
#' @return A data frame of edges with their PPIs and weights.
#' @export
edge_table <- function(fit, threshold = 0.5, predictor_names = NULL) {
  nburn   <- fit$nburn
  nsample <- fit$nsample
  q       <- ncol(fit$post_gamma)
  p       <- ncol(fit$X)
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names)) predictor_names <- paste0("x", 1:p)

  G_post    <- fit$post_G[post_idx, , drop = FALSE]
  Beta_post <- fit$post_Beta[post_idx, , drop = FALSE]

  edge_ppi    <- colMeans(G_post)
  edge_weight <- colMeans(Beta_post)

  # Reconstruct pairs from upper triangle (q x q matrix)
  pairs <- which(upper.tri(matrix(0, q, q)), arr.ind = TRUE)

  df <- data.frame(
    Node1  = pairs[, 1],
    Node2  = pairs[, 2],
    PPI    = round(edge_ppi, 4),
    Weight = round(edge_weight, 4)
  )

  # Map node indices to predictor names (for non-categorical, node index = predictor index)
  if (max(pairs) <= p) {
    df$Name1 <- predictor_names[df$Node1]
    df$Name2 <- predictor_names[df$Node2]
  }

  df$Selected <- df$PPI > threshold
  df <- df[order(-df$PPI), ]
  rownames(df) <- NULL

  return(df)
}


#' GP Component Summary
#'
#' Summarizes the GP weights (phi) and precision (r) for each predictor.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param predictor_indices Which predictors to summarize (default: selected).
#' @param credible_level Credible interval level (default 0.95).
#' @param predictor_names Optional character vector of predictor names.
#' @param covariate_names Optional character vector of covariate names.
#'
#' @return A data frame with GP parameter summaries.
#' @export
gp_summary <- function(fit, predictor_indices = NULL,
                       credible_level = 0.95,
                       predictor_names = NULL,
                       covariate_names = NULL) {

  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(fit$X)
  K       <- ncol(fit$Z)
  alpha   <- 1 - credible_level
  post_idx <- (nburn + 1):(nburn + nsample)

  if (is.null(predictor_names))  predictor_names  <- paste0("x", 1:p)
  if (is.null(covariate_names))  covariate_names  <- paste0("z", 1:K)

  if (is.null(predictor_indices)) {
    ppi <- colMeans(fit$post_gamma[post_idx, 1:p, drop = FALSE])
    predictor_indices <- which(ppi > 0.5)
  }

  results <- list()
  for (j in predictor_indices) {
    # r_j summary
    r_samples <- fit$post_r[post_idx, j]
    r_row <- data.frame(
      Predictor = predictor_names[j],
      Parameter = "r_j",
      Covariate = NA_character_,
      Mean      = mean(r_samples),
      SD        = sd(r_samples),
      CI_lower  = quantile(r_samples, probs = alpha / 2),
      CI_upper  = quantile(r_samples, probs = 1 - alpha / 2)
    )
    results[[length(results) + 1]] <- r_row

    # phi_{jk} summaries
    for (k in 1:K) {
      w_samples <- fit$post_w[j, k, post_idx]
      w_row <- data.frame(
        Predictor = predictor_names[j],
        Parameter = "phi_jk",
        Covariate = covariate_names[k],
        Mean      = mean(w_samples),
        SD        = sd(w_samples),
        CI_lower  = quantile(w_samples, probs = alpha / 2),
        CI_upper  = quantile(w_samples, probs = 1 - alpha / 2)
      )
      results[[length(results) + 1]] <- w_row
    }
  }

  df <- do.call(rbind, results)
  rownames(df) <- NULL
  return(df)
}


#' Plot ROC Curve
#'
#' Plots the Receiver Operating Characteristic curve from predicted probabilities.
#'
#' @param y_true Binary response vector (0/1).
#' @param y_prob Predicted probabilities.
#' @param method_name Label for the legend (default: "VEL-BMGM").
#'
#' @return A ggplot object.
#' @export
plot_roc <- function(y_true, y_prob, method_name = "VEL-BMGM") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for this function")
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("pROC is required for this function")

  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)

  df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )

  ggplot2::ggplot(df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
    ggplot2::annotate("text", x = 0.6, y = 0.2,
                      label = paste0(method_name, " (AUC = ", auc_val, ")"),
                      size = 4, color = "steelblue") +
    ggplot2::labs(x = "False Positive Rate", y = "True Positive Rate",
                  title = "ROC Curve") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal()
}


#' Print Method for VEL-BMGM Fit Objects
#'
#' @param x A fit object from bmgm_GP().
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns x.
#' @export
print.vel_bmgm <- function(x, ...) {
  nburn   <- x$nburn
  nsample <- x$nsample
  p       <- ncol(x$X)
  K       <- ncol(x$Z)
  n       <- nrow(x$X)
  post_idx <- (nburn + 1):(nburn + nsample)

  cat("VEL-BMGM Model Fit\n")
  cat("-------------------\n")
  cat("  Observations:      ", n, "\n")
  cat("  Predictors (p):    ", p, "\n")
  cat("  Covariates (K):    ", K, "\n")
  cat("  MCMC:              ", nburn, "burn-in +", nsample, "samples\n\n")

  # Predictor selection
  ppi <- colMeans(x$post_gamma[post_idx, 1:p, drop = FALSE])
  pnames <- if (!is.null(colnames(x$X))) colnames(x$X) else paste0("X", 1:p)
  selected <- which(ppi > 0.5)

  cat("  Selected predictors (PPI > 0.5):\n")
  if (length(selected) > 0) {
    for (j in selected) {
      cat(sprintf("    %s: PPI = %.3f\n", pnames[j], ppi[j]))
    }
  } else {
    cat("    (none)\n")
  }

  # Graph edges
  if (!is.null(x$adj_G)) {
    adj <- x$adj_G
    n_edges <- sum(adj[upper.tri(adj)] != 0)
    cat(sprintf("\n  Graph edges detected: %d\n", n_edges))
  }

  cat("\n  Use summary_bmgm(fit) for detailed output.\n")
  cat("  Use pip_table(fit) for predictor/covariate selection tables.\n")
  invisible(x)
}


#' Plot Method for VEL-BMGM Fit Objects
#'
#' Dispatches to the appropriate VEL-BMGM visualization based on the type argument.
#'
#' @param x A fit object from bmgm_GP().
#' @param type Character; plot type. One of "coefficients" (default), "ppi", "graph", "trace".
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A ggplot object (or base plot for "graph").
#' @export
plot.vel_bmgm <- function(x, type = c("coefficients", "ppi", "graph", "trace"), ...) {
  type <- match.arg(type)
  switch(type,
         coefficients = plot_varying_coefficients(x, ...),
         ppi          = plot_pip_heatmap(x, ...),
         graph        = plot_graph(x, ...),
         trace        = plot_trace(x, ...))
}
