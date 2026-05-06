# Bayesian Mixed Graphical Model with Varying Coefficient Logistic Regression

Fits a Bayesian mixed graphical model for logistic regression, allowing
predictor effects to vary with covariates via Gaussian Process priors.
Supports mixed-type predictors (continuous, discrete, zero-inflated, and
categorical) and estimates a sparse conditional independence graph with
spike-and-slab priors.

## Usage

``` r
bmgm_GP(
  X,
  Y,
  Z,
  type_y = "b_new",
  type,
  nburn = 1000,
  nsample = 1000,
  theta_priors,
  v_0 = 0.05,
  v_1 = 1,
  pi_beta,
  seed,
  context_spec = T,
  bfdr = 0.05,
  cont = FALSE,
  a = -2.75,
  b = 0.5,
  a_0 = 1,
  b_0 = 1,
  a_r = 1,
  b_r = 100,
  a_tau_w = 1,
  b_tau_w = 1,
  a_tau = 1,
  b_tau = 1,
  alpha_jk = 0.5,
  lengthscale = 0.5,
  sigma_kernel,
  tune = 100,
  kernel = "sqexp",
  ...
)
```

## Arguments

- X:

  Matrix of predictors (n x p).

- Y:

  Binary response vector (length n).

- Z:

  Matrix of covariates for varying effects (n x K).

- type_y:

  Character; response type (default: "b_new" for logistic regression).

- type:

  Vector indicating type of each predictor in X (e.g. "c" for
  continuous, "d" for discrete, "z" for zero-inflated, "m" for
  categorical).

- nburn:

  Number of burn-in MCMC iterations.

- nsample:

  Number of MCMC samples to retain.

- theta_priors:

  List of priors for each node (predictor).

- v_0, v_1:

  Variances for spike-and-slab priors.

- pi_beta:

  Bernoulli prior for variable selection (default: 2/(p-1)).

- seed:

  Optional random seed.

- context_spec:

  Logical; whether to use context-specific graph for categorical
  variables.

- bfdr:

  Bayesian FDR for edge selection.

- cont:

  Logical; apply F(X) transformation to continuous variables.

- a, b, a_0, b_0, a_r, b_r, a_tau_w, b_tau_w, a_tau, b_tau, alpha_jk,
  lengthscale, sigma_kernel, tune:

  Hyperparameters for priors and MCMC (see details).

- kernel:

  Character; kernel type for GP (default: "sqexp").

- ...:

  Additional arguments (currently ignored).

## Value

A list containing posterior samples, estimated adjacency matrices, and
fitted model parameters.
