# Selection Metrics

Computes predictor selection accuracy given true predictor indices.

## Usage

``` r
selection_metrics(fit, true_predictors, threshold = 0.5)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- true_predictors:

  Integer vector of true predictor indices.

- threshold:

  PPI threshold for selection (default 0.5).

## Value

A named vector of selection metrics.
