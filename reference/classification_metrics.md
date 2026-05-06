# Classification Metrics

Computes standard binary classification metrics from model predictions.

## Usage

``` r
classification_metrics(y_true, y_pred_prob, threshold = 0.5)
```

## Arguments

- y_true:

  True binary labels.

- y_pred_prob:

  Predicted probabilities (from predict_jbmgm).

- threshold:

  Classification threshold (default 0.5).

## Value

A named vector of metrics.
