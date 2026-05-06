# Plot Method for VEL-BMGM Fit Objects

Dispatches to the appropriate VEL-BMGM visualization based on the type
argument.

## Usage

``` r
# S3 method for class 'vel_bmgm'
plot(x, type = c("coefficients", "ppi", "graph", "trace"), ...)
```

## Arguments

- x:

  A fit object from bmgm_GP().

- type:

  Character; plot type. One of "coefficients" (default), "ppi", "graph",
  "trace".

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A ggplot object (or base plot for "graph").
