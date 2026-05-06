# Compute Cj Variance and XCXt Matrix

Computes covariance and XCXt for one variable's GP component.

## Usage

``` r
compute_Cj_var(Kernels, Xj, gamma_tilde, w, r)
```

## Arguments

- Kernels:

  Array of kernel matrices (n x n x K).

- Xj:

  Vector of values for variable j.

- gamma_tilde:

  Vector of inclusion indicators.

- w:

  Vector of kernel weights.

- r:

  Precision parameter.

## Value

List with C_j (GP covariance) and XCXt.
