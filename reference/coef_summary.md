# Coefficient Summary Table

Summarizes varying coefficients beta_j(Z) — posterior means and credible
intervals.

## Usage

``` r
coef_summary(
  fit,
  predictor_indices = NULL,
  credible_level = 0.95,
  predictor_names = NULL
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

## Value

A data frame with pointwise summaries for each observation.
