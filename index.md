## VEL.BMGM ![](reference/figures/logo.png)

**Varying-Effects Logistic Regression & Bayesian Mixed Graphical
Structure**

This R package implements Bayesian varying-effects logistic regression
for binary outcomes, with mixed-type predictors (continuous, discrete,
zero-inflated, and categorical). Predictor effects can vary nonlinearly
with external covariates, modeled via Gaussian Process priors. The
package also estimates a sparse conditional independence graph among
predictors, enabling structure discovery and variable selection.

## Key Features

- **Flexible logistic regression**: Models complex predictor effects
  using Gaussian Processes.
- **Mixed data support**: Handles continuous, discrete, zero-inflated,
  and categorical predictors.
- **Spike-and-slab priors**: For variable and covariate selection, and
  for graphical model edges.
- **Efficient MCMC**: Uses Polya-Gamma augmentation for fast posterior
  sampling.
- **Graph estimation**: Jointly estimates a sparse conditional
  independence graph among predictors.

## Installation

`VEL.BMGM` depends on the `BMGM` package, which must be installed first:

``` r

devtools::install_github("mauroflorez/BMGM")
devtools::install_github("mauroflorez/VEL.BMGM")
```

## Basic Usage

``` r

library(VEL.BMGM)
set.seed(1)

n <- 300
X <- matrix(rnorm(n * 4), n, 4)         # 4 continuous predictors
Z <- matrix(runif(n, -1, 1), n, 1)      # 1 continuous covariate

# Only X[, 1] is a true predictor; its effect varies sinusoidally with Z[, 1]
eta <- X[, 1] * 2 * sin(pi * Z[, 1])
Y   <- rbinom(n, 1, plogis(eta))

fit <- bmgm_GP(X, Y, Z, type = rep("c", 4),
               nburn = 5000, nsample = 5000, seed = 1)

# Posterior inclusion probabilities (X[, 1] should be near 1, others near 0)
colMeans(fit$post_gamma[-(1:5000), ])

# Visualize the estimated varying coefficient beta_1(Z)
plot(fit)
```

The Gaussian Process kernel can be selected via the `kernel` argument:
`"sqexp"` (the default) is the standard squared exponential kernel and
gives smooth posterior sample paths, while `"exp"` produces less smooth
paths and can be useful for capturing sharper transitions. The
`lengthscale` argument controls the smoothness of the GP; smaller values
allow more local variation.

For a more comprehensive demonstration with multiple varying-coefficient
shapes (sine, quadratic, linear, sigmoid) recovered simultaneously, see
the [recovery of varying coefficients
vignette](https://mauroflorez.github.io/VEL.BMGM/articles/recovery-of-varying-effects.md).

## Authors

Mauro Florez

This package is under active development. Please see the paper for full
model and algorithm details.
