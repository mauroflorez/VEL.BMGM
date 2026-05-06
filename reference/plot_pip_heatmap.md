# Plot PPI Heatmap

Heatmap of two-level variable selection results. Rows are predictors,
columns are covariates, color encodes covariate-level PPI.

## Usage

``` r
plot_pip_heatmap(
  fit,
  threshold = 0.5,
  predictor_names = NULL,
  covariate_names = NULL
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold for annotation (default 0.5).

- predictor_names:

  Optional character vector of predictor names.

- covariate_names:

  Optional character vector of covariate names.

## Value

A ggplot heatmap object.
