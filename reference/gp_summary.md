# GP Component Summary

Summarizes the GP weights (phi) and precision (r) for each predictor.

## Usage

``` r
gp_summary(
  fit,
  predictor_indices = NULL,
  credible_level = 0.95,
  predictor_names = NULL,
  covariate_names = NULL
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- predictor_indices:

  Which predictors to summarize (default: selected).

- credible_level:

  Credible interval level (default 0.95).

- predictor_names:

  Optional character vector of predictor names.

- covariate_names:

  Optional character vector of covariate names.

## Value

A data frame with GP parameter summaries.
