# Compute Sigma_Y for Prediction

Computes the posterior covariance matrix for Y in the Bayesian model.

## Usage

``` r
compute_Sigma_Y(XCXt_all, gamma, tau, omega)
```

## Arguments

- XCXt_all:

  Array of XCXt matrices for all variables (n x n x p).

- gamma:

  Inclusion indicators for variables.

- tau:

  Variance or precision for Y.

- omega:

  Observation-specific weights (default is 1).

## Value

Covariance matrix for Y.
