#' Squared Exponential Kernel
#'
#' Computes the squared exponential (RBF) kernel matrix for a given covariate vector.
#'
#' @param Z Covariate vector.
#' @param lengthscale Lengthscale parameter (default: 0.5).
#' @param sigma_f Output scale parameter (default: 1).
#' @param kernel Character; kernel type. "sqexp" for squared exponential (default), "exp" for exponential.
#'
#' @return Kernel matrix.
#' @export
kernel_se <- function(Z, lengthscale = 0.5, sigma_f = 1, kernel = "sqexp") {
  if (kernel == "sqexp") {
    dists <- as.matrix(dist(Z))^2
    sigma_f^2 * exp(-0.5 * dists / lengthscale^2)
  } else {
    dists <- as.matrix(dist(Z))
    sigma_f^2 * exp(-dists / lengthscale^2)
  }
}

#' Compute Cj Variance and XCXt Matrix
#'
#' Computes covariance and XCXt for one variable's GP component.
#'
#' @param Kernels Array of kernel matrices (n x n x K).
#' @param Xj Vector of values for variable j.
#' @param gamma_tilde Vector of inclusion indicators.
#' @param w Vector of kernel weights.
#' @param r Precision parameter.
#'
#' @return List with C_j (GP covariance) and XCXt.
#' @export
compute_Cj_var <- function(Kernels, Xj, gamma_tilde, w, r){
  n <- dim(Kernels)[1]
  K <- dim(Kernels)[3]
  C_j <- Reduce('+', lapply(1:K, function(k) {gamma_tilde[k]*(w[k]^2)*Kernels[,,k]})) + (1/r)*diag(n)
  XCXt <- diag(Xj) %*% C_j %*% diag(Xj)
  return(list(C_j, XCXt))
}

#' Compute Sigma_Y for Prediction
#'
#' Computes the posterior covariance matrix for Y in the Bayesian model.
#'
#' @param XCXt_all Array of XCXt matrices for all variables (n x n x p).
#' @param gamma Inclusion indicators for variables.
#' @param tau Variance or precision for Y.
#' @param omega Observation-specific weights (default is 1).
#'
#' @return Covariance matrix for Y.
#' @export
compute_Sigma_Y <- function(XCXt_all, gamma, tau, omega){
  n <- length(omega)
  included_vars <- which(gamma == 1)
  Sigma_Y <- diag(1/omega) + (1/tau)*matrix(1, nrow = n, ncol = n)
  if(length(included_vars) > 0){
    Sigma_Y <- Sigma_Y + Reduce('+', lapply(included_vars, function(j) XCXt_all[,,j]))
  }
  return(Sigma_Y)
}
