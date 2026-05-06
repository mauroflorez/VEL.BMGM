# Posterior Predictive Check

Generates posterior predictive samples and computes Bayesian p-values.

## Usage

``` r
posterior_predictive_check(fit, n_rep = 500, summary_fn = mean)
```

## Arguments

- fit:

  Model fit object from bmgm_GP().

- n_rep:

  Number of posterior predictive replications (default 500).

- summary_fn:

  Summary statistic function (default: mean).

## Value

A list with replicated stats, observed stat, and Bayesian p-value.
