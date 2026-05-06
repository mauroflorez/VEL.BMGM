# Package index

## Main fitting function

- [`bmgm_GP()`](https://mauroflorez.github.io/VEL.BMGM/reference/bmgm_GP.md)
  : Bayesian Mixed Graphical Model with Varying Coefficient Logistic
  Regression

## Prediction

- [`predict_jbmgm()`](https://mauroflorez.github.io/VEL.BMGM/reference/predict_jbmgm.md)
  : Predict with VEL-BMGM Model

## Post-fit summaries

- [`summary_bmgm()`](https://mauroflorez.github.io/VEL.BMGM/reference/summary_bmgm.md)
  : Summary of VEL-BMGM Model Fit
- [`pip_table()`](https://mauroflorez.github.io/VEL.BMGM/reference/pip_table.md)
  : Posterior Inclusion Probability Table
- [`coef_summary()`](https://mauroflorez.github.io/VEL.BMGM/reference/coef_summary.md)
  : Coefficient Summary Table
- [`edge_table()`](https://mauroflorez.github.io/VEL.BMGM/reference/edge_table.md)
  : Edge Posterior Inclusion Probability Table
- [`gp_summary()`](https://mauroflorez.github.io/VEL.BMGM/reference/gp_summary.md)
  : GP Component Summary

## Plots

- [`plot_varying_coefficients()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_varying_coefficients.md)
  : Plot Varying Coefficients
- [`plot_pip_heatmap()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_pip_heatmap.md)
  : Plot PPI Heatmap
- [`plot_pip_bar()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_pip_bar.md)
  : Plot PPI Bar Chart
- [`plot_graph()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_graph.md)
  : Plot Predictor Graph
- [`plot_trace()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_trace.md)
  : MCMC Traceplots
- [`plot_roc()`](https://mauroflorez.github.io/VEL.BMGM/reference/plot_roc.md)
  : Plot ROC Curve

## Metrics

- [`classification_metrics()`](https://mauroflorez.github.io/VEL.BMGM/reference/classification_metrics.md)
  : Classification Metrics
- [`selection_metrics()`](https://mauroflorez.github.io/VEL.BMGM/reference/selection_metrics.md)
  : Selection Metrics
- [`posterior_predictive_check()`](https://mauroflorez.github.io/VEL.BMGM/reference/posterior_predictive_check.md)
  : Posterior Predictive Check

## Internals (kernel and MCMC helpers)

- [`kernel_se()`](https://mauroflorez.github.io/VEL.BMGM/reference/kernel_se.md)
  : Squared Exponential Kernel
- [`compute_Cj_var()`](https://mauroflorez.github.io/VEL.BMGM/reference/compute_Cj_var.md)
  : Compute Cj Variance and XCXt Matrix
- [`compute_Sigma_Y()`](https://mauroflorez.github.io/VEL.BMGM/reference/compute_Sigma_Y.md)
  : Compute Sigma_Y for Prediction
- [`update_beta_j()`](https://mauroflorez.github.io/VEL.BMGM/reference/update_beta_j.md)
  : Update Beta_j for Bayesian Logistic Model
- [`update_beta_0()`](https://mauroflorez.github.io/VEL.BMGM/reference/update_beta_0.md)
  : Update Intercept Beta_0 for Bayesian Logistic Model
- [`update_tau_w()`](https://mauroflorez.github.io/VEL.BMGM/reference/update_tau_w.md)
  : Update Tau_w Hyperparameter for GP Priors

## S3 methods

- [`print(`*`<vel_bmgm>`*`)`](https://mauroflorez.github.io/VEL.BMGM/reference/print.vel_bmgm.md)
  : Print Method for VEL-BMGM Fit Objects
- [`plot(`*`<vel_bmgm>`*`)`](https://mauroflorez.github.io/VEL.BMGM/reference/plot.vel_bmgm.md)
  : Plot Method for VEL-BMGM Fit Objects
