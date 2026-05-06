# Plot Predictor Graph

Visualize the estimated predictor graph as a network diagram.

## Usage

``` r
plot_graph(
  fit,
  threshold = 0.5,
  predictor_names = NULL,
  node_color = NULL,
  layout = "circle"
)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- threshold:

  PPI threshold for edge inclusion (default 0.5).

- predictor_names:

  Optional character vector of predictor names.

- node_color:

  Optional color vector for nodes (length p).

- layout:

  igraph layout name (default "circle").

## Value

An igraph plot (invisible igraph object).
