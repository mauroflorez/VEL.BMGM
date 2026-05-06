# Summary of VEL-BMGM Model Fit

Prints a concise summary of the fitted model including predictor
selection, covariate selection, intercept, and graph structure.

## Usage

``` r
summary_bmgm(
  fit,
  threshold = 0.5,
  credible_level = 0.95,
  predictor_names = NULL,
  covariate_names = NULL
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold for selection (default 0.5).

- credible_level:

  Credible interval level (default 0.95).

- predictor_names:

  Optional character vector of predictor names.

- covariate_names:

  Optional character vector of covariate names.

## Value

Invisibly returns a list with summary components.
