# Update Intercept Beta_0 for Bayesian Logistic Model

Update Intercept Beta_0 for Bayesian Logistic Model

## Usage

``` r
update_beta_0(X, z_y, beta_j_current, gamma, omega, tau)
```

## Arguments

- X:

  Design matrix of predictors.

- z_y:

  Adjusted response vector.

- beta_j_current:

  Current beta_j matrix.

- gamma:

  Inclusion indicators for predictors.

- omega:

  Observation weights.

- tau:

  Precision for prior on beta_0.

## Value

Updated beta_0 scalar.
