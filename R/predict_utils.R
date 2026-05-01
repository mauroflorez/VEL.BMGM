#' Predict with VEL-BMGM Model
#'
#' Generates predictions for new test data from a fitted VEL-BMGM model.
#' Uses GP conditional prediction to extrapolate varying coefficients
#' from training covariate values to test covariate values.
#'
#' @param fit Model fit object from bmgm_GP().
#' @param X_test Test predictor matrix (n_test x p).
#' @param Z_test Test covariate matrix (n_test x K).
#' @param Z_train Training covariate matrix (n_train x K). If NULL, uses fit$Z.
#' @param mcmc_samples Number of posterior samples to use (default 1000).
#' @param threshold PPI threshold for covariate selection (default 0.5).
#'
#' @return A list with:
#'   \item{prob_mean}{Vector of mean predicted probabilities (length n_test).}
#'   \item{prob_samples}{Matrix (mcmc_samples x n_test) of predicted probabilities per MCMC iteration.}
#'   \item{class}{Predicted binary class (threshold at 0.5).}
#' @export
predict_jbmgm <- function(fit, X_test, Z_test, Z_train = NULL,
                          mcmc_samples = 1000, threshold = 0.5) {

  nburn   <- fit$nburn
  nsample <- fit$nsample
  p       <- ncol(X_test)
  n_test  <- nrow(X_test)
  K       <- ncol(Z_test)

  if (is.null(Z_train)) Z_train <- fit$Z
  n_train <- nrow(Z_train)

  # Normalize covariates to [0,1] using TRAINING range
  Z_all   <- rbind(Z_train, Z_test)
  Z_norm  <- apply(Z_all, 2, function(x) scales::rescale(x, to = c(0, 1)))
  Z_train_norm <- Z_norm[1:n_train, , drop = FALSE]
  Z_test_norm  <- Z_norm[(n_train + 1):(n_train + n_test), , drop = FALSE]

  # Build kernel matrices
  Kernels_train      <- array(0, dim = c(n_train, n_train, K))
  Kernels_test_train <- array(0, dim = c(n_test, n_train, K))

  lengthscale <- 0.5  # matches the fitting function
  sigma_f     <- 1

  for (k in 1:K) {
    Kernels_train[,,k] <- kernel_se(Z_train_norm[, k],
                                    lengthscale = lengthscale,
                                    sigma_f = sigma_f)

    dists_sq <- outer(Z_test_norm[, k], Z_train_norm[, k],
                      function(a, b) (a - b)^2)
    Kernels_test_train[,,k] <- sigma_f^2 * exp(-0.5 * dists_sq / lengthscale^2)
  }

  # GP conditional prediction: beta_j_star = C(Z*, Z) C(Z, Z)^{-1} beta_j_hat
  calc_Cj_pred <- function(w, r) {
    CZZstar <- Reduce('+', lapply(1:K, function(k)
      w[k]^2 * Kernels_test_train[,,k]))
    CZZ <- Reduce('+', lapply(1:K, function(k)
      w[k]^2 * Kernels_train[,,k])) + (1/r) * diag(n_train)
    CZZ_inv <- solve(CZZ)
    return(CZZstar %*% CZZ_inv)
  }

  # Compute posterior mean beta_j (training) and covariate PPIs
  post_idx <- (nburn + 1):(nburn + nsample)
  gamma_post <- fit$post_gamma[post_idx, 1:p, drop = FALSE]

  beta_j_est <- matrix(0, nrow = n_train, ncol = p)
  PPI_Z      <- matrix(0, nrow = p, ncol = K)

  for (j in 1:p) {
    included_iters <- which(gamma_post[, j] == 1)
    if (length(included_iters) > 0) {
      beta_samples <- fit$beta_j[1:n_train, j, post_idx[included_iters]]
      if (is.matrix(beta_samples)) {
        beta_j_est[, j] <- rowMeans(beta_samples)
      } else {
        beta_j_est[, j] <- beta_samples
      }
    }

    for (k in 1:K) {
      PPI_Z[j, k] <- mean(fit$post_gamma_tilde[j, k, post_idx])
    }
  }

  # Sample posterior iterations
  n_avail    <- min(mcmc_samples, nsample)
  sample_idx <- sample(1:nsample, n_avail)

  beta_j_star <- array(0, dim = c(n_test, p, n_avail))
  eta_i       <- matrix(0, nrow = n_avail, ncol = n_test)

  for (m in seq_along(sample_idx)) {
    s <- sample_idx[m] + nburn
    for (j in 1:p) {
      if (fit$post_gamma[s, j] == 1) {
        ind <- as.numeric(PPI_Z[j, ] > threshold)
        w_t <- fit$post_w[j, , s] * ind
        r_t <- fit$post_r[s, j]
        beta_j_star[, j, m] <- calc_Cj_pred(w_t, r_t) %*% beta_j_est[, j]
      }
    }

    beta_0_m <- fit$post_beta0[s]
    eta_i[m, ] <- beta_0_m + rowSums(X_test * beta_j_star[,, m])
  }

  prob_samples <- plogis(eta_i)  # logistic transform (base R, no dependency)
  prob_mean    <- colMeans(prob_samples)
  y_class      <- as.integer(prob_mean > 0.5)

  list(
    prob_mean    = prob_mean,
    prob_samples = prob_samples,
    class        = y_class
  )
}
