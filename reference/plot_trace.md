# MCMC Traceplots

Produces traceplots for key model parameters.

## Usage

``` r
plot_trace(
  fit,
  parameter = "beta0",
  predictor_indices = NULL,
  include_burnin = TRUE
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- parameter:

  One of "beta0", "gamma", "tauw", "r".

- predictor_indices:

  For gamma or r: which predictors (default: first 4).

- include_burnin:

  Whether to include burn-in iterations (default TRUE).

## Value

A ggplot object.
