# Plot Varying Coefficients

Plots estimated varying-coefficient functions beta_j(Z) with credible
bands.

## Usage

``` r
plot_varying_coefficients(
  fit,
  predictor_indices = NULL,
  covariate_index = 1,
  credible_level = 0.95,
  predictor_names = NULL,
  covariate_names = NULL,
  true_functions = NULL
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- predictor_indices:

  Which predictors to plot (default: all selected with PPI \> 0.5).

- covariate_index:

  Which covariate dimension to use as x-axis.

- credible_level:

  Credible interval level (default 0.95).

- predictor_names:

  Optional character vector of predictor names.

- covariate_names:

  Optional character vector of covariate names.

- true_functions:

  Optional named list of functions for overlay (e.g., list("x1" =
  function(z) 2\*z^2+0.5)).

## Value

A ggplot object (one facet per predictor).
