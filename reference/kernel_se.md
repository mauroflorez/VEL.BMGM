# Squared Exponential Kernel

Computes the squared exponential (RBF) kernel matrix for a given
covariate vector.

## Usage

``` r
kernel_se(Z, lengthscale = 0.5, sigma_f = 1, kernel = "sqexp")
```

## Arguments

- Z:

  Covariate vector.

- lengthscale:

  Lengthscale parameter (default: 0.5).

- sigma_f:

  Output scale parameter (default: 1).

- kernel:

  Character; kernel type. "sqexp" for squared exponential (default),
  "exp" for exponential.

## Value

Kernel matrix.
