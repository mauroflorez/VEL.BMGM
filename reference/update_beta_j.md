# Update Beta_j for Bayesian Logistic Model

Updates the beta_j coefficient vector for a single predictor using
Bayesian updates.

## Usage

``` r
update_beta_j(X, j, z_y, Cj_var, gamma, omega_pg, beta_0, beta_j_current)
```

## Arguments

- X:

  Design matrix of predictors.

- j:

  Index of predictor to update.

- z_y:

  Adjusted response vector.

- Cj_var:

  List or array of GP covariance matrices.

- gamma:

  Inclusion indicators for predictors.

- omega_pg:

  Polya-Gamma weights.

- beta_0:

  Intercept.

- beta_j_current:

  Current beta_j matrix.

## Value

Updated beta_j vector for predictor j.
