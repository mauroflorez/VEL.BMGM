# Posterior Inclusion Probability Table

Produces a publication-ready data frame of PPIs at both selection
levels.

## Usage

``` r
pip_table(fit, threshold = 0.5, predictor_names = NULL, covariate_names = NULL)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold for selection (default 0.5).

- predictor_names:

  Optional character vector of predictor names.

- covariate_names:

  Optional character vector of covariate names.

## Value

A data frame with predictor and covariate PPIs.
