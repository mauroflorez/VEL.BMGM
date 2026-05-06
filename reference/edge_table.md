# Edge Posterior Inclusion Probability Table

Produces a table of edge PPIs for the predictor graph.

## Usage

``` r
edge_table(fit, threshold = 0.5, predictor_names = NULL)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold for edge inclusion (default 0.5).

- predictor_names:

  Optional character vector of predictor names.

## Value

A data frame of edges with their PPIs and weights.
