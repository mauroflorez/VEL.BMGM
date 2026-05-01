## VEL.BMGM <img src="man/figures/logo.png" align="right" width="100"/>

**Varying-Effects Logistic Regression & Bayesian Mixed Graphical Structure**

This R package implements Bayesian varying-effects logistic regression for binary outcomes, with mixed-type predictors (continuous, discrete, zero-inflated, and categorical). Predictor effects can vary nonlinearly with external covariates, modeled via Gaussian Process priors. The package also estimates a sparse conditional independence graph among predictors, enabling structure discovery and variable selection.

## Key Features

-   **Flexible logistic regression**: Models complex predictor effects using Gaussian Processes.
-   **Mixed data support**: Handles continuous, discrete, zero-inflated, and categorical predictors.
-   **Spike-and-slab priors**: For variable and covariate selection, and for graphical model edges.
-   **Efficient MCMC**: Uses Polya-Gamma augmentation for fast posterior sampling.
-   **Graph estimation**: Jointly estimates a sparse conditional independence graph among predictors.

## Installation

`VEL.BMGM` depends on the `BMGM` package, which must be installed first:

``` r
devtools::install_github("mauroflorez/BMGM")
devtools::install_github("mauroflorez/VEL-BMGM")
```

## Basic Usage

``` r
library(VEL.BMGM)
fit <- bmgm_GP(X, Y, Z, type_y = 'b_new', type = type_vector)
```

The Gaussian Process kernel can be selected via the `kernel` argument
(`"sqexp"`, the default, or `"exp"` for rougher sample paths that better
capture sharp non-monotone effects).

## Notes for users upgrading from 0.1.0

The squared-exponential kernel formula was changed to the standard form
(with the 1/2 factor in the exponent). At the same nominal `lengthscale`,
the new version produces smoother varying-coefficient curves than the old
one. If you previously tuned `lengthscale`, expect to re-tune.

## Authors

Mauro Florez

This package is under active development. Please see the paper for full model and algorithm details.
