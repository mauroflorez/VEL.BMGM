# Plot ROC Curve

Plots the Receiver Operating Characteristic curve from predicted
probabilities.

## Usage

``` r
plot_roc(y_true, y_prob, method_name = "VEL-BMGM")
```

## Arguments

- y_true:

  Binary response vector (0/1).

- y_prob:

  Predicted probabilities.

- method_name:

  Label for the legend (default: "VEL-BMGM").

## Value

A ggplot object.
