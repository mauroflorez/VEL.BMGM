#' Bayesian Mixed Graphical Model with Varying Coefficient Logistic Regression
#'
#' Fits a Bayesian mixed graphical model for logistic regression, allowing predictor effects to vary with covariates via Gaussian Process priors. Supports mixed-type predictors (continuous, discrete, zero-inflated, and categorical) and estimates a sparse conditional independence graph with spike-and-slab priors.
#'
#' @param X Matrix of predictors (n x p).
#' @param Y Binary response vector (length n).
#' @param Z Matrix of covariates for varying effects (n x K).
#' @param type_y Character; response type (default: "b_new" for logistic regression).
#' @param type Vector indicating type of each predictor in X (e.g. "c" for continuous, "d" for discrete, "z" for zero-inflated, "m" for categorical).
#' @param nburn Number of burn-in MCMC iterations.
#' @param nsample Number of MCMC samples to retain.
#' @param theta_priors List of priors for each node (predictor).
#' @param v_0,v_1 Variances for spike-and-slab priors.
#' @param pi_beta Bernoulli prior for variable selection (default: 2/(p-1)).
#' @param seed Optional random seed.
#' @param context_spec Logical; whether to use context-specific graph for categorical variables.
#' @param bfdr Bayesian FDR for edge selection.
#' @param cont Logical; apply F(X) transformation to continuous variables.
#' @param a,b,a_0,b_0,a_r,b_r,a_tau_w,b_tau_w,a_tau,b_tau,alpha_jk,lengthscale,sigma_kernel,tune Hyperparameters for priors and MCMC (see details).
#' @param kernel Character; kernel type for GP (default: "sqexp").
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list containing posterior samples, estimated adjacency matrices, and fitted model parameters.
#' @export

bmgm_GP <- function(X, Y, Z, type_y = 'b_new', type, nburn = 1000, nsample = 1000, theta_priors,
                            v_0 = 0.05, v_1 = 1, pi_beta, seed, context_spec = T,
                            bfdr = 0.05, cont = FALSE, a = -2.75, b = 0.5, a_0 = 1, b_0 = 1, a_r = 1, b_r = 100,
                            a_tau_w = 1, b_tau_w = 1, a_tau = 1, b_tau = 1, alpha_jk = 0.5,
                            lengthscale = 0.3, sigma_kernel, tune = 100, kernel = "sqexp",...){
  
  #X: matrix of predictors (n*p)
  #Y: response variable (n*1)
  #Z: matrix of covariates (n*K)
  #type_y: class of y (continuous(c), categorical(m))
  #type : vector indicating the type of predictor (X) ('c', 'd', 'm')
  #nburn: nburn
  #nsample: nsample
  #theta_priors: list where each element corresponds to the priors for each node (predictor)
  #v_0, v_1: variance on continuos spike-slab for graph selection
  #pi_beta: bernoulli prior on spike-slab (2/p-1) default
  #bfdr: bayesian FDR
  #cont: indicator if transformation F(X) is applied also to continuous
  #context_spec: type of graph for categorical variables
  #a, b: scalar hyperparameters for MRF prior
  #a_r, b_r  For r_j 
  #a_tau_w, b_tau_w for variance of phi in GP
  #a_tau, b_tau: for variance of response (tau)
  #alpha_jk: prior mean bernoulli for covariate selection
  
  if(!missing(seed)) set.seed(seed)
  
  if(missing(type)){
    stop("Error: type of variables is missing")
  }
  
  
  X_input <- X
  Z_input <- Z
  
  Z <- apply(Z, 2, function(x) scales::rescale(x, to = c(0, 1)))
  #Sample size n
  n <- nrow(X)
  #Number of predictors: p
  p <- ncol(X)
  #Number of covariates:
  K <- ncol(Z)
  #n
  n_i = rep(1,n)
  
  #Initial Imputation
  R <- 1 - is.na(X)*1
  r_imp <- which(rowSums(R) < p)
  means <- colMeans(X, na.rm = T)
  
  if(length(r_imp) > 0) {
    for(s in 1:p){
      impute_value <- switch(type[s],
                             "c" = means[s],
                             "d" = round(means[s]),
                             "z" = round(means[s]),
                             "m" = as.numeric(names(which.max(table(X[, s])))))
      X[, s][R[, s] == 0] <- impute_value
    }
  }
  
  #Centering of continuous
  X[, type == "c"] <- scale(X[, type == "c"], center = means[type == "c"], scale = F)
  
  #Spike and slab
  log_prior_beta <- function(Beta, pi_beta, v0, v1){
    sum(log(pi_beta*exp(dnorm(Beta, mean = 0, sd = v1, log = TRUE)) +
              (1 - pi_beta)*exp(dnorm(Beta, mean = 0, sd = v0, log = TRUE))))
  }
  
  #Split categorical variables in columns (if any):
  split <- BMGM::split_X_cat(X, type)
  X_design <- split$matrix
  categories <- split$categories
  
  #Number of nodes: q
  q <- ncol(X_design)
  #Type of each node
  type_q <- rep(type, categories)
  
  #############======= Auxiliary functions to get the names =======#############
  
  #Get the names of the nodes
  var_names <- rep(1:p, categories)
  tags <- BMGM::get_names_graph(p, q, categories, var_names)
  
  #Get the subnames for categorical nodes
  tag <- tags$tags
  #Get the names for the edges
  tag_Beta <- tags$tag_Beta
  #Get the indicators where there are no edges (among sub-categories)
  ind_noedge <- tags$indicators_noedge
  
  ####################======= Transformation ========###########################
  
  lambda <- BMGM::find_lambda(X, type)
  
  F_X <- BMGM::F_transformation(X = X_design, type = type_q, parameter = lambda, cont)
  std_err <- apply(F_X, 2, sd)
  F_scaled <- scale(F_X, center = F, scale = std_err)
  
  #Beta
  F_centered <- scale(F_scaled, center = T, scale = F)
  S <- crossprod(F_centered)
  pdxid <- diag(solve(var(F_scaled)))
  
  
  #pdxid <- rep(1,p)
  ########################===== Hyperpriors =====##############################
  
  #Theta_priors should be a list where each element corresponds to the priors
  #for each variable.
  #In the case of categorical data simply the uniform. By the default are given
  #as below, should be specified in specific scenarios.
  
  theta_priors <- list()
  for (s in 1:p) {
    theta_priors[[s]] <- switch(type[s],
                                "c" = c(mean_mu = 0.001, mean_sd = 0.001, sd_shape = 0.001, sd_rate = 0.001),
                                "d" = c(mu_shape = 0.001, mu_rate = 0.001, nu_shape = 0.001, nu_rate = 0.001),
                                "z" = c(p_shape1 = 0.001, p_shape2 = 0.001, mu_shape = 0.001, mu_rate = 0.001),
                                "m" = rep(1 / (categories[s] + 1), categories[s] + 1))
  }
  
  ##### Priors for Beta / Now as inputs
  
  if(missing(pi_beta)) pi_beta <- 2/(p-1)
  
  #hyperparameters regression
  
  if(missing(sigma_kernel)) sigma_kernel <- rep(1, K)
  #h_0 = 100
  #h_bet = max(var(X)) #within the range of variability of X  
  #b_0 = 100 #error term
  #a_0 = 1 #error term
  
  #a = -2.75 #MRF
  #b = 0.5 #MRF
  
  
  ########################===== Initial Values =====############################
  
  # Initial Values - MCMC Algorithm
  Beta <- diag(pdxid)
  G <- diag(q) - diag(q)
  theta <- list()
  gamma = rep(0, q)
  gamma_tilde = w = matrix(0, nrow = q, ncol = K)
  
  tau_w = 1
  candr = 100
  tau_y = 1
  beta_0 = 0
  tau_0 = 1
  
  n_y = rep(1, n)
  kappa_y = Y - n_y/2
  
  beta_j  = array(0, c(n, q, nburn + nsample))
  
  for(s in 1:p){
    theta[[s]] <- switch(type[s],
                         'c' = c(0,0.1),
                         'd' = c(0.1, 1),
                         'z' = c(0.5, 1),
                         'm' = rep(1/(categories[s]+1), categories[s]+1))
  }
  
  # Initial variance for steps in MCMC
  h_beta <- rep(0.1, q)
  h_theta <- list()
  scale_theta <- c()
  for(s in 1:p){
    h_theta[[s]] <- switch(type[s],
                           'c' = diag(0.1, 2),
                           'd' = diag(c(0.01, 0.01)),
                           'm' = diag(0.1, categories[s]+1),
                           'z' = diag(c(0.005, 0.5)))
    
    scale_theta[s] <- switch(type[s],
                             'c' = 2.4^2/2,
                             'd' = 2.4^2/2,
                             'z' = 2.4^2/2,
                             'm' = 2.4^2/(categories[s]+1))
  }
  
  #New parameters - GP:
  
  Kernels_Z = array(0, dim = c(n, n, K))
  Cj_var <- array(0, dim = c(n, n, q))
  XCXt_var <- array(0, dim = c(n, n, q))
  #Cj_beta <- array(0, dim = c(n, n, q))
  
  for(k in 1:K){
    Kernels_Z[,,k] <- kernel_se(Z[,k], lengthscale = lengthscale,
                                sigma_f = sigma_kernel[k], kernel = kernel)
  }
  
  #Since we are fixing sigma_kernel and tau_w, we need to calculate the
  #Kernels just one time.
  
  lastgammak_all = matrix(0, nrow = q, ncol = K)
  lastw_all = matrix(0, nrow = q, ncol = K)
  lastr_all = rep(100, q)
  
  for(j in 1:q){
    Cj_var[,,j] <- diag(n)/lastr_all[j]
  }
  
  D = 1
  
  proposeWBetween = rep(0, q)
  proposeWWithin = rep(0, q)
  countWBetween = rep(0, q)
  countWWithin = rep(0, q)
  
  lambda_g = alpha_jk
  zetaProposal_all = 0.5*matrix(1, nrow = q, ncol = K)
  updateCount = rep(1, q)
  lastUpdateNum = rep(1, q)
  updateInterval = 20;
  firstUpdate = updateInterval;
  stepIncrement = 0.3
  iter_add = rep(1, q)
  m_update = 100
  
  r <- rep(1, q)
  
  n_add_prop = n_remove_prop = n_keep_prop = 0
  
  # Current values:
  beta_j_current = beta_j[,,1]
  
  ac_anterior = ac_Beta = ac_gamma = ac_theta = rep(0, q)
  n_gamma_prop =  n_keep_prop =  n_add_prop = n_remove_prop = 0
  n_gamma_accept = n_add_accept = n_remove_accept = n_keep_accept = 0
  
  ############################===== Storage =====###############################
  
  post_Beta <- matrix(nrow = nburn+nsample, ncol = q*(q-1)/2)
  post_G <- matrix(nrow = nburn+nsample, ncol = q*(q-1)/2)
  post_imputation <- list()
  post_theta <- list()
  
  #GP:
  post_omegapg <- matrix(0, nrow = nburn+nsample, ncol = n)
  post_gamma <- matrix(nrow = nburn+nsample, ncol = q)
  post_gamma_tilde = array(0, c(p, K, nburn+nsample))
  post_r = matrix(1, nrow = nburn+nsample, ncol = p)
  post_w = array(0, c(p, K, nburn+nsample))
  post_tau <- rep(0, nburn+nsample)
  post_tau0 <- rep(0, nburn+nsample)
  post_beta0 <- rep(0, nburn+nsample)
  post_tauw <- rep(0, nburn+nsample)
  gammak_temp = array(0,c(q,K,nburn+nsample));
  
  for(s in 1:p){
    post_theta[[s]] <- switch(type[s],
                              'c' = matrix(nrow = nburn+nsample, ncol = 2),
                              'd' = matrix(nrow = nburn+nsample, ncol = 2),
                              'z' = matrix(nrow = nburn+nsample, ncol = 2),
                              'm' = matrix(nrow = nburn+nsample, ncol = categories[s]+1))
  }
  
  colnames(post_G) <- tag_Beta[upper.tri(tag_Beta)]
  colnames(post_Beta) <- tag_Beta[upper.tri(tag_Beta)]
  
  ############################ Normalization Constants ############################
  
  log_norm_constant <- function(s, parameters){
    domain <- 0:100
    se <- std_err[s]
    un_llk_z <- function(x, theta)
      exp(theta[2]*(log(theta[1])*x - lfactorial(x)) - theta[3]*BMGM::F_transformation(x, type = "d", lambda)/se)
    
    log(apply(parameters, 1, function(x) sum(un_llk_z(domain, theta = x))))
  }
  
  log_norm_constant_Z <- function(s, parameters){
    domain <- 0:100
    se <- std_err[s]
    un_llk <- function(x, theta){
      (theta[1]*(x == 0) + (1 - theta[1])*dpois(x, lambda = theta[2]))*
        exp(-theta[3]*BMGM::F_transformation(x, type = "d", lambda)/se)
    }
    log(apply(parameters, 1, function(x) sum(un_llk(domain, theta = x))))
  }
  
  #######################======= MCMC Algorithm =======##########################
  
  #pb <- txtProgressBar(min = 0, max = nburn + nsample, style = 3)
  
  for(m in 1:(nburn + nsample)){
    ######################## 1st Block - Updates Theta ###########################
    for(s in 1:p){
      theta_s <- theta[[s]]
      h_theta_s <- h_theta[[s]]
      
      #Proposal
      switch(type[s],
             'c' = {
               #sampling mu
               mu0 = theta_priors[[s]][1]
               tau0 = theta_priors[[s]][2]
               
               mu = means[s]#mean(X[,s])
               tao = theta_s[2]
               
               var_post <- tau0 + n*tao
               mean_post <- (tau0*mu0 + n*tao*mu)/var_post
               
               #adding norm. constant
               C_s = sum(c(F_scaled[,-s]%*%Beta[s,-s]))
               
               mean_post <- mean_post - (1/var_post)*C_s
               theta_s[1] <- rnorm(1, mean = mean_post, sd = 1/sqrt(var_post))
               
               #sampling tau. In this case we need MCMC
               
               tau <- theta_s[2]
               mu <- theta_s[1]
               tau_proposal <- abs(rnorm(1, mean = tau, h_theta_s))
               
               a0 = theta_priors[[s]][3]
               b0 = theta_priors[[s]][4]
               
               shape_post = a0 + n/2
               rate_post = b0 + 0.5*sum((X[,s] - mu)^2)
               
               C_s_2 = sum(c(F_scaled[,-s]%*%Beta[s,-s])^2)
               
               ar <- dgamma(tau_proposal, shape = shape_post,rate = rate_post, log = T) -
                 dgamma(tau, shape = shape_post, rate = rate_post, log = T) +
                 C_s_2/2*(1/tau - 1/tau_proposal)
               
               accept <- min(1, exp(ar))
               
               if(stats::runif(1) < accept){
                 theta_s[2] <- tau_proposal
                 ac_theta[s] <- ac_theta[s] + 1
               }
               
             },
             'd' = {
               theta_star <- abs(MASS::mvrnorm(n = 1, mu = c(theta_s[1], theta_s[2]),
                                               Sigma = h_theta_s))
               
               C_s = c(F_scaled[,-s]%*%Beta[s,-s])
               
               param_s <- cbind('mu' = rep(theta_s[1], n), 'nu' = rep(theta_s[2], n), 'edge' = C_s)
               param_star <- cbind('mu' = rep(theta_star[1], n), 'nu' = rep(theta_star[2], n), 'edge' = C_s)
               
               log_Z_s <- sum(log_norm_constant(s, param_s))
               log_Z_star <- sum(log_norm_constant(s, param_star))
               
               log_density_d <- sum(theta_s[2]*(log(theta_s[1])*X[,s] - lfactorial(X[,s])))
               log_density_d_star <- sum(theta_star[2]*(log(theta_star[1])*X[,s] - lfactorial(X[,s])))
               
               theta_priors_s <- theta_priors[[s]]
               
               log_prior_d <- dgamma(theta_s[1], shape = theta_priors_s[1], rate = theta_priors_s[2], log = TRUE) +
                 dgamma(theta_s[2], shape = theta_priors_s[3], rate = theta_priors_s[4], log = TRUE)
               
               log_prior_d_star <- dgamma(theta_star[1], shape = theta_priors_s[1], rate = theta_priors_s[2], log = TRUE) +
                 dgamma(theta_star[2], shape = theta_priors_s[3], rate = theta_priors_s[4], log = TRUE)
               
               #acceptance:
               
               log_ar <- log_density_d_star - log_density_d +
                 log_prior_d_star - log_prior_d +
                 log_Z_s - log_Z_star
               
               accept <- min(1, exp(log_ar))
               
               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }
               
             },
             'z' = {
               theta_star <- c(mnormt::rmtruncnorm(n = 1, mean = theta_s[1], varcov = h_theta_s[1,1], lower = 0, upper = 1),
                               abs(rnorm(n = 1, mean = theta_s[2], sd = h_theta_s[2,2])))
               
               C_s <- c(F_scaled[,-s]%*%Beta[s,-s])
               
               prior_values <- theta_priors[[s]]
               alpha0 <- prior_values[1] #assc with p
               beta0 <- prior_values[2]
               
               a0 <- prior_values[3] #with mu
               b0 <- prior_values[4]
               
               m_Z = sum(X[,s] == 0)
               
               post_s1 <- m_Z + alpha0
               post_s2 <- n - m_Z + beta0
               
               post_alpha <- a0 + sum(X[,s])
               post_beta <- n - m_Z + b0
               
               #Normalizing constant
               
               param_s <- cbind('pi' = rep(theta_s[1], n), 'mu' = rep(theta_s[2], n), 'edge' = C_s)
               param_star <- cbind('pi' = rep(theta_star[1], n), 'mu' = rep(theta_star[2], n), 'edge' = C_s)
               
               log_Z_s <- sum(log_norm_constant_Z(s, param_s))
               log_Z_star <- sum(log_norm_constant_Z(s, param_star))
               
               log_ar <- dbeta(theta_star[1], post_s1, post_s2, log = T) - dbeta(theta_s[1], post_s1, post_s2, log = T) +
                 dgamma(theta_star[2], post_alpha, post_beta, log = T) - dgamma(theta_s[2], post_alpha, post_beta, log = T) +
                 log_Z_s - log_Z_star
               
               accept <- min(1, exp(log_ar))
               
               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }
             },
             'm' = {
               theta_star <- abs(MASS::mvrnorm(n = 1, mu = theta_s, Sigma = h_theta_s))
               theta_star <- theta_star/sum(theta_star)
               
               cat = as.numeric(factor(X[,s]))
               #edge-potentials
               cols = which(var_names == s)
               se <- std_err[cols]
               C_s <- F_scaled[,-cols]%*%Beta[-cols,cols]
               
               un_llk <- cbind(rep(exp(log(theta_s[1])),n),
                               exp(apply(C_s, 1, function(x) log(theta_s[-1]) - x/se)))
               norm_consts <- rowSums(un_llk)
               
               un_llk_star <-  cbind(rep(exp(log(theta_star[1])),n),
                                     (exp(apply(C_s, 1, function(x) log(theta_star[-1]) - x/se))))
               norm_consts_star <- rowSums(un_llk_star)
               
               llk  <- un_llk[cbind(1:nrow(C_s), cat)]/norm_consts
               llk_star  <- un_llk_star[cbind(1:nrow(C_s), cat)]/norm_consts
               
               log_dif_llk <- sum(log(llk_star) - log(llk))
               
               prior_s <- theta_priors[[s]]
               #priors
               log_dif_priors <- log(gtools::ddirichlet(theta_star, prior_s)) -
                 log(gtools::ddirichlet(theta_s, prior_s))
               
               log_ar <- log_dif_llk + log_dif_priors
               
               accept <- min(1, exp(log_ar))
               
               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }
             })
      
      theta[[s]] <- theta_s
      post_theta[[s]][m,] <- theta[[s]]
    }
    
    ########################### 2nd Block - Update Beta ##########################
    # Arkaprava & Dunson (2020)
    
    #update = "beta"
    Beta_star = Beta
    for(l in 1:q){
      theta_s <- theta[[var_names[l]]]
      mean <- S[l,-l] #mean
      Omega <- Beta[-l, -l]
      diag(Omega) <- pmax(pdxid[-l], 1e-6)
      Omegai <- eigen(Omega)
      OmegatempiU <- t(Omegai$vectors)/sqrt(abs(Omegai$values))
      
      #Update column l
      Omega_inv <- crossprod(OmegatempiU)
      
      Ci <- eigen((S[l,l] + 1)*Omega_inv + diag(ifelse(G[l,-l] == 0, 1/v_0, 1/v_1)))
      CiU <- t(Ci$vectors)/sqrt(abs(Ci$values))
      C_inv <- crossprod(CiU)
      
      #Proposal
      mean_proposal <- -C_inv%*%mean
      var_proposal <- C_inv
      
      Beta_proposal <- MASS::mvrnorm(1, mean_proposal, var_proposal)
      Beta_proposal[which(is.na(Beta_proposal))] <- 0
      
      #Adjust the update wrt the acceptance rate
      k_2 <- min(1, as.numeric(h_beta[l]/sqrt(crossprod(Beta_proposal - Beta[l,-l]))))
      Beta_proposal <- Beta[l,-l] + k_2*(Beta_proposal - Beta[l, -l])
      
      Beta_star[l, -l] <- Beta_proposal
      Beta_star[-l, l] <- Beta_proposal
      Beta_star <- Beta_star*ind_noedge
      
      C_s = c(F_scaled[,-l]%*%Beta[l,-l])
      C_star = c(F_scaled[,-l]%*%Beta_proposal)
      
      if(type[var_names[l]] == "c"){
        mu = theta_s[1]
        tau = theta_s[2]
        
        mean = mu - C_s/tau
        mean_star = mu - C_star/tau
        
        log_dif_llk <- sum(dnorm(X[,l], mean = mean_star, sd = sqrt(1/tau), log = T)) -
          sum(dnorm(X[,l], mean = mean, sd = sqrt(1/tau), log = T))
      } else {
        log_dif_llk <- sum(F_scaled[,l]*C_star) - sum(F_scaled[,l]*C_s)
      }
      
      log_dif_priors <- log_prior_beta(Beta_proposal, pi_beta, v_0, v_1) -
        log_prior_beta(Beta[l,-l], pi_beta, v_0, v_1)
      
      log_dif_prop <- mvtnorm::dmvnorm(Beta[-l, l], mean_proposal, var_proposal, log = T) -
        mvtnorm::dmvnorm(Beta_proposal, mean_proposal, var_proposal, log = T)
      
      switch(type[var_names[l]],
             'd' = {
               param_s <- cbind('mu' = rep(theta_s[1], n),
                                'nu' = rep(theta_s[2], n),
                                'edge' = C_s)
               param_star <- cbind('mu' = rep(theta_s[1], n),
                                   'nu' = rep(theta_s[2], n),
                                   'edge' = C_star)
               
               log_Z_s <- sum(log_norm_constant(s, param_s))
               log_Z_star <- sum(log_norm_constant(s, param_star))
               log_dif_norm <- log_Z_s - log_Z_star
             },
             'c' = {
               mu = theta_s[1]
               tau = theta_s[2]
               log_Z_s <- sum(mu*tau*C_s - C_s^2/(2*tau))
               log_Z_star <- sum(mu*tau*C_star - C_star^2/(2*tau))
               log_dif_norm = log_Z_s - log_Z_star
             },
             'z' = {
               param_s <- cbind('p' = rep(theta_s[1], n),
                                'mu' = rep(theta_s[2], n),
                                'edge' = C_s)
               param_star <- cbind('p' = rep(theta_s[1], n),
                                   'mu' = rep(theta_s[2], n),
                                   'edge' = C_star)
               
               log_Z_s <- sum(log_norm_constant_Z(s, param_s))
               log_Z_star <- sum(log_norm_constant_Z(s, param_star))
               log_dif_norm <- log_Z_s - log_Z_star
             },
             'm' = {
               Beta_star <- Beta
               Beta_star[l, -l] = Beta_star[-l,l] = Beta_proposal
               cat = as.numeric(factor(X[,s]))
               #edge-potentials
               cols = which(var_names == s)
               se <- std_err[cols]
               C_s <- F_scaled[,-cols]%*%Beta[-cols,cols]
               C_star <- F_scaled[,-cols]%*%Beta_star[-cols,cols]
               
               un_llk <- cbind(rep(exp(log(theta_s[1])),n),
                               exp(apply(C_s, 1, function(x) log(theta_s[-1]) - x/se)))
               log_Z_s <- sum(log(rowSums(un_llk)))
               
               un_llk_star <-  cbind(rep(exp(log(theta_s[1])),n),
                                     exp(apply(C_star, 1, function(x) log(theta_s[-1]) - x/se)))
               log_Z_star <- sum(log(rowSums(un_llk_star)))
               log_dif_norm <- log_Z_s - log_Z_star
             })
      
      log_ar <- log_dif_llk + log_dif_priors + log_dif_norm
      
      accept <- min(1, exp(log_ar))
      
      if(stats::runif(1) < accept){
        Beta = Beta_star
        ac_Beta[l] <- ac_Beta[l] + 1
      }
    }
    
    ######################== 3rd Block - Update G and pi ===###################
    G_new <- Beta[upper.tri(Beta)]
    slab <- pi_beta*dnorm(G_new, 0, v_1)
    spike <- (1-pi_beta)*dnorm(G_new, 0, v_0)
    G_new <- slab/(slab+spike)
    nan_ind <- is.nan(G_new)
    G_new[nan_ind] <- 0
    
    G_vector <- rbinom(q*(q-1)/2, size = 1, prob = G_new)
    
    G <- matrix(0, ncol = q, nrow = q)
    G[upper.tri(Beta)] <- G_vector
    G <- G + t(G)
    
    post_Beta[m,] <- Beta[upper.tri(Beta)]
    post_G[m,] <- G_vector
    
    ########################=== 4th Block - Regression ===######################
    
    #If the response is continuous:
    
    if(type_y == 'c'){
      ## 1. Update tau:
      #Proposal
      tau_try <- rgamma(n = 1, 2, tau_y/2)
      #Evaluate
      prior_ratio <- -(a_tau + 1)*log(tau_try) - b_tau/tau_try + ((a_tau + 1)*log(tau_y) + b_tau/tau_y)
      
      var_try = llk_Y(X, Z, tau_try, gamma, w, sigma_kernel, r)
      var_y = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel, r)
      
      llk_try <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_try[[1]], log = T) 
      llk <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_y[[1]], log = T) 
      
      qlastnew = (tune - 1)*log(tau_try) - tune*tau_try/tau_y 
      
      qnewlast = (tune - 1)*log(tau_y) - tune*tau_y/tau_try 
      
      alpha = prior_ratio + qnewlast - qlastnew + llk_try - llk
      
      if(log(runif(1)) < alpha){
        tau_y = tau_try
      }
      
      ## 2. Update Gamma:
      
      change_index = sample(1:q, 1)
      gamma_prop = gamma
      
      if(gamma[change_index] == 0){
        #Add new variable
        gamma_prop[change_index] = 1
        n_add_prop <- n_add_prop + 1
      } else {
        keep_or_drop = rbinom(1, 1, 0.5)
        #Keep or eliminate
        gamma_prop[change_index] = keep_or_drop
        
        if(keep_or_drop == 1) n_keep_prop <- n_keep_prop + 1
        else n_remove_prop <- n_remove_prop + 1
      }
      
      j = change_index
      k = sample(1:K, 1)
      
      last_r = post_r[max(m - 1, 1),j]
      last_w = post_w[j,k,m]
      
      
      if(gamma_prop[j] == 1){ #we added or kept
        
        if(gamma_tilde[j,k] == 1){
          gamma_tilde_prop = 0
          w_prop = 0
        }
        else {
          gamma_tilde_prop = 1
          w_prop = rnorm(1, mean = 0, sd = tau_w)
        }
        
        #If keep we use the previous values for sampling, if no, sample from priors
        if(gamma[change_index] == 0){
          #added
          r_prop = rgamma(1, a_r, 1/b_r)
          sigma_kernel_prop = rgamma(1, a_kernel, 1/b_kernel)
          
          qlastnew = qnewlast = 0
          
          indicator = 1
        } else {
          #kept
          r_prop = rgamma(1, tune, last_r/tune)
          sigma_kernel_prop = rgamma(1, tune, last_sigma_kernel/tune)
          
          qlastnew = (tune - 1)*log(r_prop) - tune*r_prop/last_r +
            (tune - 1)*log(sigma_kernel_prop) - tune*sigma_kernel_prop/last_sigma_kernel
          
          qnewlast = (tune - 1)*log(last_r) - tune*last_r/r_prop +
            (tune - 1)*log(last_sigma_kernel) - tune*last_sigma_kernel/sigma_kernel_prop
          indicator = 2
        }
      } else {
        #removed
        gamma_tilde_prop = 0
        w_prop = 0
        r_prop = rgamma(1, a_r, 1/b_r)
        sigma_kernel_prop = rgamma(1, a_kernel, 1/b_kernel)
        
        qlastnew = qnewlast = 0
        indicator = 0
      }
      
      r_try = r
      r_try[j] = r_prop
      gamma_tilde_try = gamma_tilde
      gamma_tilde_try[j,k] = gamma_tilde_prop
      gamma_try = gamma_prop
      w_try = w
      w_try[j,k] = w_prop
      sigma_kernel_try = sigma_kernel
      sigma_kernel_try[k] <- sigma_kernel_prop
      
      var_try = llk_Y(X, Z, tau_y, gamma_try, w_try, sigma_kernel_try, r_try)
      var_y = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel, r)
      
      llk_try <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_try[[1]], log = T) 
      llk <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_y[[1]], log = T) 
      
      prior_try <- dgamma(r_prop, a_r, 1/b_r, log = T) + 
        dgamma(sigma_kernel_prop, a_kernel, 1/b_kernel, log = T)  
      
      prior <- dgamma(last_r, a_r, 1/b_r, log = T) + 
        dgamma(last_sigma_kernel, a_kernel, 1/b_kernel, log = T)
      
      mrf_gamma_try = a*sum(gamma_try) + b*t(gamma_try)%*%G%*%gamma_try 
      
      mrf_gamma = a*sum(gamma) + b*t(gamma)%*%G%*%gamma
      #G should be 0 in the diagonal
      
      ar <- llk_try - llk + prior_try - prior + qnewlast - qlastnew + mrf_gamma_try - mrf_gamma
      
      if(log(runif(1)) < ar){
        #Accept
        gamma <- gamma_prop
        gamma_tilde <- gamma_tilde_try
        w <- w_try
        r[j] <- r_prop
        sigma_kernel <- sigma_kernel_try
        gamma_tilde <- gamma_tilde_try
        
        if(indicator == 0) n_remove_accept = n_remove_accept + 1
        else if(indicator == 1) n_add_accept = n_add_accept + 1
        else if(indicator == 2) n_keep_accept = n_keep_accept + 1
      }
      
      #Level 2:
      
      if(gamma[j] == 1){
        #if we accepted
        for(k in 1:K){
          #level 2 - between
          
          if(gamma_tilde[j,k] == 1){
            gamma_tilde_prop = 0
            w_prop = 0
          } else {
            gamma_tilde_prop = 1
            w_prop = rnorm(1, mean = 0, sd = tau_w)
          }
          w_try = w
          w_try[j,k] = w_prop
          
          var_try = llk_Y(X, Z, tau_y, gamma, w_try, sigma_kernel, r)
          var_y = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel, r)
          
          llk_try <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_try[[1]], log = T) 
          llk <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_y[[1]], log = T) 
          
          prior_try <- (gamma_tilde_prop*dnorm(w_prop, mean = 0, sd = tau_w) + 
                          (1 - gamma_tilde_prop)*(w_prop == 0))*dbinom(gamma_tilde_prop, 1, 0.5)
          
          prior <- gamma_tilde[j,k]*dnorm(w[j,k], mean = 0, sd = tau_w) +
            (1 - gamma_tilde[j,k])*(w[j,k] == 0)*dbinom(gamma_tilde[j,k], 1, 0.5)
          
          ar <- llk_try - llk + log(prior_try) - log(prior)
          
          if(log(runif(1)) < ar){
            gamma_tilde[j,k] = gamma_tilde_prop
            w[j,k] = w_prop
          }
          
          #level 2 within:
          
          if(gamma_tilde_prop == 1){
            w_prop = rnorm(1, mean = 0, sd = tau_w)
            
            w_try = w
            w_try[j,k] = w_prop
            
            var_try = llk_Y(X, Z, tau_y, gamma, w_try, sigma_kernel, r)
            var_y = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel, r)
            
            llk_try <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_try[[1]], log = T) 
            llk <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_y[[1]], log = T) 
            
            prior_try <- dnorm(w_prop, mean = 0, sd = tau_w)*dbinom(1, 1, 0.5)
            prior <- dnorm(w[j,k], mean = 0, sd = tau_w)*dbinom(1, 1, 0.5)
            
            ar <- llk_try - llk + log(prior_try) - log(prior)
            
            if(log(runif(1)) < ar){
              gamma_tilde[j,k] = 1
              w[j,k] = w_prop
            }
          }
          
          #update r, sigma_kernel
          last_sigma_kernel = sigma_kernel[k]
          last_r = r[j]
          
          r_prop = rgamma(1, tune, last_r/tune)
          sigma_kernel_prop = rgamma(1, tune, last_sigma_kernel/tune)
          
          qlastnew = (tune - 1)*log(r_prop) - tune*r_prop/last_r +
            (tune - 1)*log(sigma_kernel_prop) - tune*sigma_kernel_prop/last_sigma_kernel
          
          qnewlast = (tune - 1)*log(last_r) - tune*last_r/r_prop +
            (tune - 1)*log(last_sigma_kernel) - tune*last_sigma_kernel/sigma_kernel_prop
          
          r_try = r
          r_try[j] = r_prop
          sigma_kernel_try = sigma_kernel
          sigma_kernel_try[k] <- sigma_kernel_prop
          
          var_try = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel_try, r_try)
          var_y = llk_Y(X, Z, tau_y, gamma, w, sigma_kernel, r)
          
          llk_try <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_try[[1]], log = T) 
          llk <- mvnfast::dmvn(Y, mu = rep(0,n), sigma = var_y[[1]], log = T) 
          
          prior_try <- dgamma(r_prop, a_r, 1/b_r, log = T) + 
            dgamma(sigma_kernel_prop, a_kernel, 1/b_kernel, log = T)  
          
          prior <- dgamma(last_r, a_r, 1/b_r, log = T) + 
            dgamma(last_sigma_kernel, a_kernel, 1/b_kernel, log = T)
          
          ar <- llk_try - llk + prior_try - prior + qnewlast - qlastnew
          
          if(log(runif(1)) < ar){
            #Accept
            r[j] <- r_prop
            sigma_kernel[k] <- sigma_kernel_prop
          }
        }
      }
      
      #Draw Beta_j(Z)
      if(gamma[j] == 1){
        var_beta = llk_Y(X_design, Z, tau_y, gamma, w, sigma_kernel, r)
        
        sigma_tau = 1/(tau_y^2)*diag(n)
        sigma_j = solve(solve(var_beta[[2]][,,j]) + t(diag(X[,j]))%*%solve(sigma_tau)%*%diag(X[,j]))
        
        sum_h = 0
        for(h in (1:q)[-j]){
          sum_h <- sum_h + gamma[h]*X[,h]*beta_j_current[,h]
        }
        mu_j = sigma_j%*%solve(sigma_tau)%*%diag(X[,j])%*%(Y - sum_h)
        
        beta_j_current[,j] = c(mvnfast::rmvn(1, mu = mu_j, sigma = sigma_j))
      } else{
        beta_j_current[,j] = rep(0,n)
      }
      
      post_gamma[m, ] <- gamma
      post_gamma_tilde[,,m] <- gamma_tilde  
      post_r[m,] <- r
      post_w[,,m] <- w
      post_sigma_kernel[m,] <- sigma_kernel
      beta_j[,,m] <- beta_j_current
      post_tau[m] <- tau_y
      
    } else if(type_y == 'b'){
      #Binomial case:
      
      change_index = sample(1:q, 1)
      j = change_index
      
      #n_gamma_prop = n_gamma_prop + 1
      
      gamma_prop = gamma
      gamma_prop[change_index] = abs(gamma[change_index] - 1)
      
      if(gamma[change_index] == 0){
        #Add new variable
        gamma_prop[change_index] = 1
        n_add_prop <- n_add_prop + 1
      } else {
        keep_or_drop = rbinom(1, 1, 0.5)
        #Keep or eliminate
        gamma_prop[change_index] = keep_or_drop
        
        if(keep_or_drop == 1) n_keep_prop <- n_keep_prop + 1
        else n_remove_prop <- n_remove_prop + 1
      }
      
      # First update omega_pg
      
      included_vars <- which(gamma == 1)
      psi_i <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
      
      omega_pg <- BayesLogit::rpg(num = n, h = n_i, z = psi_i)
      z_y <- kappa_y/omega_pg
      
      lastr = lastr_all[j]
      lastw = lastw_all[j,]
      lastgammak = lastgammak_all[j,]
      
      zetaProposal = zetaProposal_all[j,]
      
      var_y = compute_Sigma_Y(XCXt_all = XCXt_var, gamma = gamma, tau = tau_0, omega = omega_pg)
      
      llk <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_y, log = T)
      
      k = sample(1:K, 1)
      
      candgammak = lastgammak
      candw = lastw
      
      if(gamma_prop[j] == 1){ #we added or kept
        #candgammak = lastgammak
        #candw = lastw
        component = rbinom(1, 1, prob = lambda_g)
        candgammak[k] = (1-component)*rbinom(1,1, zetaProposal[k]) + 
          component*rbinom(1,1, 0.5)
        
        if(candgammak[k] == 1){
          #proposewbetween
          candw[k] = rnorm(1, mean = 0, sd = 1/sqrt(tau_w))
        } else {
          candw[k] = 0
        }
        
        if(gamma[change_index] == 0){
          indicator = 1
          candr = rgamma(1, shape = a_r, scale = b_r)
          
          qlastnew = qnewlast = 0
          
        } else {
          indicator = 2
          candr = rgamma(n = 1, shape = tune, scale = lastr/tune)
          
          qlastnew = (tune - 1)*log(candr) - tune*candr/lastr 
          qnewlast = (tune - 1)*log(lastr) - tune*lastr/candr 
        } 
      } else {
        #candgammak[k] = 0
        #candw[k] = 0
        candgammak = rep(0, K)
        candw = rep(0, K)
        candr = rgamma(n = 1, shape = a_r, scale = b_r)
        
        qlastnew = qnewlast = 0
        indicator = 0
      } 
      
      XCXt_try <- XCXt_var
      Cj_var_try <- Cj_var
      
      update_try <- compute_Cj_var(Kernels = Kernels_Z, Xj = X[,j], gamma_tilde = candgammak, w = candw, r = candr)
      Cj_var_try[,,j] <- update_try[[1]]
      XCXt_try[,,j] <- update_try[[2]]
      
      var_try = compute_Sigma_Y(XCXt_all = XCXt_try, gamma_prop, tau = tau_0, omega = omega_pg)
      llk_try <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_try, log = T) 
      
      prior_try <- dgamma(candr, shape = a_r, scale = b_r, log = T)  
      prior <- dgamma(lastr_all[j], shape = a_r, scale = b_r, log = T) 
      
      mrf_gamma_try = a*sum(gamma_prop) + b*t(gamma_prop)%*%G%*%gamma_prop 
      mrf_gamma = a*sum(gamma) + b*t(gamma)%*%G%*%gamma
      #G should be 0 in the diagonal
      
      ar <- llk_try - llk + prior_try - prior + qnewlast - qlastnew + mrf_gamma_try - mrf_gamma
      
      if(log(runif(1)) < ar){
        #Accept
        gamma <- gamma_prop
        lastgammak_all[j,] <- candgammak
        lastw_all[j,] <- candw
        lastr_all[j] <- candr
        
        Cj_var <- Cj_var_try
        XCXt_var <- XCXt_try
        llk <- llk_try
        
        if(indicator == 0) n_remove_accept = n_remove_accept + 1
        else if(indicator == 1) n_add_accept = n_add_accept + 1
        else if(indicator == 2) n_keep_accept = n_keep_accept + 1
      }
      
      #Level 2:
      
      if(gamma[j] == 1){
        #if we accepted
        for(k in 1:K){
          #level 2 - between
          lastgammak = lastgammak_all[j,]
          candgammak = lastgammak
          candw = lastw_all[j,]
          
          component = rbinom(1, 1, prob = lambda_g)
          candgammak[k] = (1-component)*rbinom(1,1, zetaProposal[k]) + 
            component*rbinom(1,1, 0.5)
          
          if(candgammak[k] == 1){
            
            candw[k] <- rnorm(1, mean = 0, sd = 1/sqrt(tau_w))
            
            XCXt_try <- XCXt_var
            Cj_var_try <- Cj_var
            
            update_try <- compute_Cj_var(Kernels = Kernels_Z, Xj = X[,j], gamma_tilde = candgammak, 
                                         w = candw, r = lastr_all[j])
            Cj_var_try[,,j] <- update_try[[1]]
            XCXt_try[,,j] <- update_try[[2]]
            
            var_try = compute_Sigma_Y(XCXt_try, gamma, tau = tau_0, omega = omega_pg)
            
            #var_try = llk_Y(X, Z, tau_0, gamma, w_try, sigma_kernel, lastr_all, omega = omega_pg)
            #var_y = llk_Y(X, Z, tau_0, gamma, lastw_all, sigma_kernel, lastr_all, omega = omega_pg)
            
            llk_try <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_try, log = T) 
            #llk <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_y, log = T) 
            
            prior_try <- (candgammak[k]*dnorm(candw[k], mean = 0, sd = 1/sqrt(tau_w)) + 
                            (1 - candgammak[k])*(candw[k] == 0))*dbinom(candgammak[k], 1, 0.5)
            
            prior <- (lastgammak_all[j,k]*dnorm(lastw_all[j,k], mean = 0, sd = 1/sqrt(tau_w)) +
                        (1 - lastgammak_all[j,k])*(lastw_all[j,k] == 0))*dbinom(lastgammak_all[j,k], 1, 0.5)
            
            ar <- llk_try - llk + log(prior_try) - log(prior)
            
            if(log(runif(1)) < ar){
              lastgammak_all[j,] <- candgammak
              lastw_all[j,] = candw
              
              Cj_var <- Cj_var_try
              XCXt_var <- XCXt_try
              llk <- llk_try
            }
            
          } else if(lastgammak[k] != 0){
            
            candw[k] <- 0
            
            XCXt_try <- XCXt_var
            Cj_var_try <- Cj_var
            
            update_try <- compute_Cj_var(Kernels = Kernels_Z, Xj = X[,j], gamma_tilde = candgammak, 
                                         w = candw, r = lastr_all[j])
            
            Cj_var_try[,,j] <- update_try[[1]]
            XCXt_try[,,j] <- update_try[[2]]
            
            var_try = compute_Sigma_Y(XCXt_try, gamma, tau = tau_0, omega = omega_pg)
            
            llk_try <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_try, log = T) 
            
            
            prior_try <- (candgammak[k]*dnorm(candw[k], mean = 0, sd = 1/sqrt(tau_w)) + 
                            (1 - candgammak[k])*(candw[k] == 0))*dbinom(candgammak[k], 1, 0.5)
            
            prior <- (lastgammak_all[j,k]*dnorm(lastw_all[j,k], mean = 0, sd = 1/sqrt(tau_w)) +
                        (1 - lastgammak_all[j,k])*(lastw_all[j,k] == 0))*dbinom(lastgammak_all[j,k], 1, 0.5)
            
            ar <- llk_try - llk + log(prior_try) - log(prior)
            
            if(log(runif(1)) < ar){
              lastgammak_all[j,] <- candgammak
              lastw_all[j,] = candw
              
              Cj_var <- Cj_var_try
              XCXt_var <- XCXt_try
              llk <- llk_try
            }
          }
          
          lastgammak = lastgammak_all[j,]
          #level 2 within:
          
          if(lastgammak[k] == 1){
            for(d in 1:D){
              
              candw[k] = rnorm(1, mean = 0, sd = 1/sqrt(tau_w))
              
              XCXt_try <- XCXt_var
              Cj_var_try <- Cj_var
              
              update_try <- compute_Cj_var(Kernels = Kernels_Z, Xj = X[,j], gamma_tilde = lastgammak, 
                                           w = candw, r = lastr_all[j])
              
              Cj_var_try[,,j] <- update_try[[1]]
              XCXt_try[,,j] <- update_try[[2]]
              
              var_try = compute_Sigma_Y(XCXt_try, gamma, tau = tau_0, omega = omega_pg)
              
              llk_try <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_try, log = T) 
              
              prior_try <- dnorm(candw[k], mean = 0, sd = 1/sqrt(tau_w))*dbinom(1, 1, 0.5)
              prior <- dnorm(lastw_all[j,k], mean = 0, sd = 1/sqrt(tau_w))*dbinom(1, 1, 0.5)
              
              ar <- llk_try - llk + log(prior_try) - log(prior)
              
              if(log(runif(1)) < ar){
                lastw_all[j,] = candw
                
                Cj_var <- Cj_var_try
                XCXt_var <- XCXt_try
                llk <- llk_try
              }
            }
          }
        }
        
        #update r, sigma_kernel
        last_r = lastr_all[j]
        
        candr = rgamma(1, shape = tune, scale = last_r/tune)
        qlastnew = (tune - 1)*log(candr) - tune*candr/last_r 
        qnewlast = (tune - 1)*log(last_r) - tune*last_r/candr 
        
        XCXt_try <- XCXt_var
        Cj_var_try <- Cj_var
        
        update_try <- compute_Cj_var(Kernels = Kernels_Z, Xj = X[,j], gamma_tilde = lastgammak_all[j,], 
                                     w = lastw_all[j,], r = candr)
        
        Cj_var_try[,,j] <- update_try[[1]]
        XCXt_try[,,j] <- update_try[[2]]
        
        var_try = compute_Sigma_Y(XCXt_try, gamma, tau = tau_0, omega = omega_pg)
        
        llk_try <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_try, log = T) 
        #llk <- mvnfast::dmvn(z_y, mu = rep(0,n), sigma = var_y, log = T) 
        
        prior_try <- dgamma(candr, shape = a_r, scale = b_r, log = T) 
        
        prior <- dgamma(last_r, shape = a_r, scale = b_r, log = T) 
        
        ar <- llk_try - llk + prior_try - prior + qnewlast - qlastnew
        
        if(log(runif(1)) < ar){
          #Accept
          lastr_all[j] <- candr
          
          Cj_var <- Cj_var_try
          XCXt_var <- XCXt_try
          llk <- llk_try
        }
      }
      
      gammak_temp[,,iter_add[j]] = lastgammak_all
      testInd = (updateCount[j] == 1)
      updateNum = updateInterval*updateCount[j]
      
      if(iter_add[j] == updateNum*(1-testInd) + firstUpdate*testInd){
        updateCount[j] = updateCount[j] + 1
        sumGammak = rowSums(gammak_temp[j,,lastUpdateNum[j]:iter_add[j], drop = F], dims = 2)
        iterations = iter_add[j] - lastUpdateNum[j] + 1
        ergAvgZeta = sumGammak/iterations
        step = stepIncrement/updateCount[j]
        zetaProposal = zetaProposal + step*(ergAvgZeta - zetaProposal)
        lastUpdateNum[j] = updateNum
      }
      
      iter_add[j] = iter_add[j] + 1
      zetaProposal_all[j,] = zetaProposal
      
      ###########################################################################
      
      #Draw Beta_j(Z)
      
      beta_j_current[,j] <- update_beta_j(X, j, z_y, Cj_var, gamma, omega_pg, beta_0, beta_j_current)
      
      beta_0 <- update_beta_0(X, z_y, beta_j_current, gamma, omega_pg, tau_0)
      
      ###########################################################################
      
      #update tau_w:
      
      tau_w <- update_tau_w(lastgammak_all, lastw_all, a_tau_w, b_tau_w)
      
      
      ### CHECK TAU0
      tau_0 <- rgamma(1, shape = a_0 + 0.5, rate = b_0 + beta_0^2/2)
      
      post_tauw[m] <- tau_w
      post_gamma[m, ] <- gamma
      post_gamma_tilde[,,m] <- lastgammak_all  
      post_r[m,] <- lastr_all
      post_w[,,m] <- lastw_all
      beta_j[,,m] <- beta_j_current
      post_tau0[m] <- tau_0
      post_beta0[m] <- beta_0
    } else if(type_y == "b_new"){
      
      #new algorithm
      
      included_vars <- which(gamma == 1)
      eta_current <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
      
      omega_pg <- BayesLogit::rpg(num = n, h = n_i, z = eta_current)
      z_y <- kappa_y/omega_pg
      
      j = sample(1:p, 1)
      
      ### update gamma:
      
      gamma_prop = gamma
      gamma_prop[j] = 1 - gamma_prop[j]
      
      # if(gamma[j] == 0){
      #   #added new
      #   gamma_prop[j] = 1
      #   n_add_prop <- n_add_prop + 1
      # } else {
      #   keep_or_drop = rbinom(1,1, 0.5)
      #   gamma_prop[j] = keep_or_drop
      # 
      #   if(keep_or_drop == 1) n_keep_prop <- n_keep_prop + 1
      #   else n_remove_prop <- n_remove_prop
      # }
      gamma_tilde_prop = gamma_tilde
      phi_prop = lastw_all
      
      zetaProposal = zetaProposal_all[j,]
      
      k <- sample(1:K, 1)
      lastr = lastr_all[j]
      
      if(gamma_prop[j] == 1){
        
        component = rbinom(1, 1, prob = lambda_g)
        gamma_tilde_prop[j,k] = (1-component)*rbinom(1,1, zetaProposal[k]) + 
          component*rbinom(1,1, 0.5)
        
        if(gamma_tilde_prop[j,k] == 1){
          phi_prop[j,k] = rnorm(1, mean = 0, sd = sqrt(1/tau_w)) 
        } else {
          phi_prop[j,k] = 0
        }
        
        #Build Cj
        Cj_var_try = Cj_var
        Cj_var_try[,,j] = Reduce('+', lapply(1:K, function(k) {
          gamma_tilde_prop[j,k]*(phi_prop[j,k]^2)*Kernels_Z[,,k]})) + 
          (1/lastr)*diag(n)
        
        #sample beta_j
        
        included_vars <- which(gamma_prop == 1 & (1:p != j))
        eta_y <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
        
        epsilon_j <- z_y - eta_y
        
        V_omega_inv <- chol(diag(X[,j]^2*omega_pg) + solve(Cj_var_try[,,j]))
        V_omega <- chol2inv(V_omega_inv)
        
        m_omega <- V_omega%*%(X[,j]*omega_pg*epsilon_j)
        
        beta_j_try = beta_j_current
        beta_j_try[,j] = mvnfast::rmvn(1, mu = m_omega, sigma = V_omega)
        #beta_j_try[,j] = mvnfast::rmvn(1, mu = rep(0,n), sigma =  Cj_var_try[,,j])
        
        included_vars <- which(gamma == 1 & 1:p != j)
        eta_minus_j <- beta_0 + rowSums(beta_j_try[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
        
        llk_0 <- mvnfast::dmvn(z_y, eta_minus_j, diag(1/omega_pg), log = TRUE)
        margllk_1 <- mvnfast::dmvn(z_y, eta_minus_j, diag(X[,j]) %*% Cj_var_try[,,j] %*% t(diag(X[,j])) + diag(1/omega_pg), log = TRUE)
        
        mrf_gamma = a * sum(gamma) + b * as.numeric(t(gamma) %*% G %*% gamma)
        mrf_gamma_try = a * sum(gamma_prop) + b * as.numeric(t(gamma_prop) %*% G %*% gamma_prop)
        
        log_prior_phi <- ifelse(gamma_tilde_prop[j,k]==1, dnorm(phi_prop[j,k],0,sqrt(1/tau_w),log=T)+log(alpha_jk), log(1-alpha_jk))
        
        log_alpha <- (margllk_1 - llk_0) + 
          (mrf_gamma_try - mrf_gamma) +
          log_prior_phi
        
      } else{
        gamma_tilde_prop = gamma_tilde
        phi_prop = lastw_all
        
        gamma_tilde_prop[j,] = rep(0,K)
        phi_prop[j,] = rep(0,K)
        
        # included_vars <- which(gamma_prop == 1 & (1:p != j))
        # eta_y <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
        # 
        Cj_var_try = Cj_var
        Cj_var_try[,,j] = (1/lastr)*diag(n)
        # 
        # V_omega_inv <- chol(diag(X[,j]^2*omega_pg) + solve(Cj_var_try[,,j]))
        # V_omega <- chol2inv(V_omega_inv)
        # 
        # m_omega <- V_omega%*%(X[,j]*omega_pg*epsilon_j)
        beta_j_try = beta_j_current
        beta_j_try[,j] = rep(0,n)
        
        
        # prior_phi <- sum(sapply(1:K, function(k) {
        #   if (gamma_tilde[j,k] == 1) dnorm(lastw_all[j,k], 0, sqrt(1/tau_w), log = TRUE) + log(alpha_jk)
        #   else log(1 - alpha_jk)
        # }))
        
        mrf_gamma = a * sum(gamma) + b * as.numeric(t(gamma) %*% G %*% gamma)
        mrf_gamma_try = a * sum(gamma_prop) + b * as.numeric(t(gamma_prop) %*% G %*% gamma_prop)
        
        included_vars <- which(gamma == 1 & 1:p != j)
        eta_minus_j <- beta_0 + rowSums(beta_j_try[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
        
        llk_0 <- mvnfast::dmvn(z_y, eta_minus_j, diag(1/omega_pg), log = TRUE)
        margllk_1 <- mvnfast::dmvn(z_y, eta_minus_j, diag(X[,j]) %*% Cj_var[,,j] %*% t(diag(X[,j])) + diag(1/omega_pg), log = TRUE)
        
        log_prior_phi <- ifelse(gamma_tilde_prop[j,k]==1, dnorm(phi_prop[j,k],0,sqrt(1/tau_w),log=T)+log(alpha_jk), log(1-alpha_jk))
        
        log_alpha <- (llk_0 - margllk_1) + 
          (mrf_gamma_try - mrf_gamma) -
          log_prior_phi
      }
      
      
      if (log(runif(1)) < log_alpha) {
        # Accept: update all relevant parameters!
        gamma <- gamma_prop
        gamma_tilde <- gamma_tilde_prop
        lastw_all <- phi_prop
        beta_j_current <- beta_j_try
        
        Cj_var = Cj_var_try
      }
      
      if(gamma[j] == 1) {
        beta_j_vec <- beta_j_current[, j]
        for (k in 1:K) {
          # Propose to toggle gamma_tilde[j,k]
          gamma_tilde_new <- gamma_tilde[j, ]
          phi_new         <- lastw_all[j, ]
          
          component = rbinom(1, 1, prob = lambda_g)
          gamma_tilde_new[k] = (1-component)*rbinom(1,1, zetaProposal[k]) + 
            component*rbinom(1,1, 0.5)
          
          if (gamma_tilde_new[k] == 1) {
            phi_new[k] <- rnorm(1, 0, sqrt(1/tau_w))
          } else {
            phi_new[k] <- 0
          }
          
          # Build new C_j (rank-1 update)
          Cj_var_try = Cj_var
          Cj_var_try[,,j] <- Reduce('+', lapply(1:K, function(l) {
            gamma_tilde_new[l] * (phi_new[l]^2) * Kernels_Z[,,l]
          })) + diag(n)/lastr_all[j]
          
          #sample new beta:
          beta_j_try <- c(update_beta_j(X, j, z_y, Cj_var_try, gamma, omega_pg, beta_0, beta_j_current))
          
          included_vars <- which(gamma == 1 & 1:p != j)
          eta_minus_j <- beta_0 + rowSums(beta_j_current[, included_vars, drop = FALSE]* X[,included_vars, drop = FALSE])
          
          margllk_0 <- mvnfast::dmvn(z_y, eta_minus_j, diag(X[,j]) %*% Cj_var[,,j] %*% t(diag(X[,j])) + diag(1/omega_pg), log = TRUE)
          margllk_1 <- mvnfast::dmvn(z_y, eta_minus_j, diag(X[,j]) %*% Cj_var_try[,,j] %*% t(diag(X[,j])) + diag(1/omega_pg), log = TRUE)

          # GP prior (multivariate normal)
          log_gp_prior_current <- mvnfast::dmvn(beta_j_vec, rep(0,n), Cj_var[,,j], log=TRUE)
          log_gp_prior_new     <- mvnfast::dmvn(beta_j_try, rep(0,n), Cj_var_try[,,j], log=TRUE)
          
          log_prior_phi_current <- sum(sapply(1:K, function(k) {
            if (gamma_tilde[j,k] == 1) dnorm(lastw_all[j,k], 0, sqrt(1/tau_w), log = TRUE) + log(alpha_jk)
            else log(1 - alpha_jk)
          }))
          
          log_prior_phi_new <- sum(sapply(1:K, function(k) {
            if (gamma_tilde_prop[j,k] == 1) dnorm(phi_new[k], 0, sqrt(1/tau_w), log = TRUE) + log(alpha_jk)
            else log(1 - alpha_jk)
          }))

          
          # MH acceptance probability
          log_alpha <- (margllk_1 - margllk_0)
            (log_gp_prior_new - log_gp_prior_current) +
            (log_prior_phi_new - log_prior_phi_current)
          
          if (log(runif(1)) < log_alpha) {
            gamma_tilde[j, k] <- gamma_tilde_new[k]
            lastw_all[j, k]   <- phi_new[k]
            
            Cj_var = Cj_var_try
          }
        }
      }
      
      for(jj in which(gamma == 1)){
        
        beta_j_vec <- beta_j_current[, jj]
        # Propose new r_j
        r0 <- lastr_all[jj]
        # Random walk on log scale: r1 = r0 * exp(N(0, s^2)), or gamma RW
        r1 <- rgamma(1, shape = tune, scale = r0 / tune)  # common choice
        
        # Build current and proposed C_j
        Cj_current <- Cj_var[,,jj]
        Cj_new <- Reduce('+', lapply(1:K, function(k)
          gamma_tilde[jj,k] * (lastw_all[jj,k]^2) * Kernels_Z[,,k]
        )) + diag(n)/r1
        
        # GP prior log-density
        log_gp_prior_current <- mvnfast::dmvn(beta_j_vec, rep(0, n), Cj_current, log=TRUE)
        log_gp_prior_new     <- mvnfast::dmvn(beta_j_vec, rep(0, n), Cj_new, log=TRUE)
        
        # Prior for r_j (gamma)
        log_prior_r0 <- dgamma(r0, shape = a_r, scale = b_r, log = TRUE)
        log_prior_r1 <- dgamma(r1, shape = a_r, scale = b_r, log = TRUE)
        
        # Proposal correction (gamma random walk, symmetric in log-scale)
        log_q_01 <- dgamma(r1, shape = tune, scale = r0/tune, log = TRUE)
        log_q_10 <- dgamma(r0, shape = tune, scale = r1/tune, log = TRUE)
        
        log_alpha <- (log_gp_prior_new - log_gp_prior_current) +
          (log_prior_r1 - log_prior_r0) +
          (log_q_10 - log_q_01)
        
        if (log(runif(1)) < log_alpha) {
          lastr_all[jj] <- r1
          # You might want to store or refresh Cj_current here for later use
          Cj_var[,,jj] <- Cj_new
        }
      }
      
    beta_j_current[,j] <- update_beta_j(X, j, z_y, Cj_var, gamma, omega_pg, beta_0, beta_j_current)
    
    zetaProposal = zetaProposal_all[j,]
    gammak_temp[,,iter_add[j]] = gamma_tilde
    testInd = (updateCount[j] == 1)
    updateNum = updateInterval*updateCount[j]
    
    if(iter_add[j] == updateNum*(1-testInd) + firstUpdate*testInd){
      updateCount[j] = updateCount[j] + 1
      sumGammak = rowSums(gammak_temp[j,,lastUpdateNum[j]:iter_add[j], drop = F], dims = 2)
      iterations = iter_add[j] - lastUpdateNum[j] + 1
      ergAvgZeta = sumGammak/iterations
      step = stepIncrement/updateCount[j]
      zetaProposal = zetaProposal + step*(ergAvgZeta - zetaProposal)
      lastUpdateNum[j] = updateNum
    }
    
    iter_add[j] = iter_add[j] + 1
    zetaProposal_all[j,] = zetaProposal

    
    ###########################################################################
    
    #Draw Beta_j(Z)
    
    #beta_j_current[,j] <- update_beta_j(X, j, z_y, Cj_var, gamma, omega_pg, beta_0, beta_j_current)
    
    beta_0 <- update_beta_0(X, z_y, beta_j_current, gamma, omega_pg, tau_0)
    
    ###########################################################################
    
    #update tau_w:
    
    tau_w <- update_tau_w(gamma_tilde, lastw_all, a_tau_w, b_tau_w)
    
    
    ### CHECK TAU0
    tau_0 <- rgamma(1, shape = a_0 + 0.5, rate = b_0 + beta_0^2/2)
    
    post_omegapg[m,] <- omega_pg
    post_tauw[m] <- tau_w
    post_gamma[m, ] <- gamma
    post_gamma_tilde[,,m] <- gamma_tilde  
    post_r[m,] <- lastr_all
    post_w[,,m] <- lastw_all
    beta_j[,,m] <- beta_j_current
    post_tau0[m] <- tau_0
    post_beta0[m] <- beta_0
  }
    
    #######################=== Adjust rates ===##############################
    
    if((m%%m_update)==0){# & m < nburn){
      ar_beta <- ac_Beta / m
      ar_theta <- ac_theta / m
      
      h_beta[ar_beta < .2] <- h_beta[ar_beta < .2]/2
      h_beta[ar_beta > .6] <- h_beta[ar_beta > .6]*2
      
      scale_theta[ar_theta < .2] <- scale_theta[ar_theta < .2]/2
      scale_theta[ar_theta > .6] <- scale_theta[ar_theta > .6]*2
      
      for(l in 1:p){
        h_theta[[l]] <- scale_theta[l]*cov(post_theta[[l]][1:m,]) +
          scale_theta[l]*0.00001*diag(length(theta[[l]]))
      }
      
      #Imputation
      for(i in r_imp){
        new_values <- BMGM::sampler_bmgm(n = 10, Beta = Beta, theta = theta, type = type,
                                   categories = categories, lambda = lambda, std = std_err,
                                   X_new = X[i,], variables = c(1 - R[i,]))
        new_values_mean <- colMeans(new_values)
        
        # Default: use mean for all
        X[i,] <- new_values_mean
        
        # Now fix "m" type (categorical)
        for (j in 1:p) {
          if (type[j] == "m" && R[i, j] == 0) {
            # Take the mode (most common category)
            X[i, j] <- as.numeric(names(which.max(table(new_values[, j]))))
          }
        }
        
        X[i, type %in% c("d", "z") & R[i,] == 0] <- round(X[i, type %in% c("d", "z") & R[i,] == 0])
      }
      post_imputation[[m/m_update]] <- X
      
      split <- BMGM::split_X_cat(X, type)
      X_design <- split$matrix
      lambda <- BMGM::find_lambda(X, type)
      
      F_X <- BMGM::F_transformation(X = X_design, type = type_q, parameter = lambda, cont)
      
      std_err <- apply(F_X, 2, sd)
      F_scaled <- scale(F_X, center = F, scale = std_err)
    }
    #pb$tick()

    print(paste0("Iteration ", m, ": |gamma| = ", sum(gamma), " / |G| = ", sum(G/2), 
                 " / Log-llk = ", round(llk_0, 1)))
  }
  
  #close(pb)
  
  ######################====== Create Adj. Matrix ======########################
  # Use only post-burn-in samples for graph estimation. Passing the full chain
  # dilutes the edge inclusion probabilities and causes the BFDR procedure to
  # drop high-confidence edges that entered the model later during burn-in.
  post_idx <- (nburn + 1):(nburn + nsample)
  ce_graph <- BMGM::context_spec_graph(q, post_Beta[post_idx, , drop = FALSE],
                                       post_G[post_idx, , drop = FALSE],
                                       tag, bfdr)
  #General Graph
  cat_graph <- BMGM::categories_graph(q, p, var_names, ce_graph$ce_esti_Z,
                                ce_graph$ce_esti_Beta, categories)
  
  #########################========= Return ==========##########################
  
  fit_our <- list(post_Beta = -post_Beta, post_theta = post_theta, post_gamma = post_gamma,
                      post_G = post_G, adj_Beta = -cat_graph$Adj_Beta, adj_G = cat_graph$Adj_Z, 
                      lambda = lambda, std = std_err, X = X_input, Z = Z_input, Y = Y, type = type, post_gamma = post_gamma, 
                      post_gamma_tilde = post_gamma_tilde, post_r = post_r, post_w = post_w, post_omegapg = post_omegapg,
                      post_tau = post_tau0, beta_j = beta_j, post_beta0 = post_beta0, post_tauw = post_tauw,
                      nburn = nburn, nsample = nsample)
  
  if(context_spec == T & any(type == "m")){
    fit_our[["adj_Beta_ce"]] <- ce_graph$ce_esti_Beta
    fit_our[["adj_Z_ce"]] <- ce_graph$ce_esti_Z
  }
  
  class(fit_our) <- "vel_bmgm"
  return(fit_our)
}
