# Plot PPI Bar Chart

Bar plot of Level-1 posterior inclusion probabilities with threshold
line.

## Usage

``` r
plot_pip_bar(fit, threshold = 0.5, predictor_names = NULL)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold (default 0.5).

- predictor_names:

  Optional character vector of predictor names.

## Value

A ggplot object.
