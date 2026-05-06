# Predict with VEL-BMGM Model

Generates predictions for new test data from a fitted VEL-BMGM model.
Uses GP conditional prediction to extrapolate varying coefficients from
training covariate values to test covariate values.

## Usage

``` r
predict_jbmgm(
  fit,
  X_test,
  Z_test,
  Z_train = NULL,
  mcmc_samples = 1000,
  threshold = 0.5
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- X_test:

  Test predictor matrix (n_test x p).

- Z_test:

  Test covariate matrix (n_test x K).

- Z_train:

  Training covariate matrix (n_train x K). If NULL, uses fit\$Z.

- mcmc_samples:

  Number of posterior samples to use (default 1000).

- threshold:

  PPI threshold for covariate selection (default 0.5).

## Value

A list with:

- prob_mean:

  Vector of mean predicted probabilities (length n_test).

- prob_samples:

  Matrix (mcmc_samples x n_test) of predicted probabilities per MCMC
  iteration.

- class:

  Predicted binary class (threshold at 0.5).
