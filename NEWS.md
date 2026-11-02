# bayesintervals 3.0.0

First CRAN release.

## New features

* `bayesian_statistical_intervals()` computes credible, prediction and
  tolerance intervals from MCMC draws of a Bayesian linear mixed model, with
  optional covariates (`covariate_cols`, `covariate_values`) and a log-normal
  back-transformation (`log_normal`).
* Two-sided tolerance intervals with asymmetric limits can be requested with
  `use_krishnamoorthy = TRUE`. `plot_krishnamoorthy()` and
  `quantile_scatterplot()` visualize the lower/upper quantile draws and the
  optimal point.
* Reducible imprecision (`ReducibleUpper`, `ReducibleLower`) is
  returned with the intervals, and `plot_one_sided_intervals()` can display
  them.
* `create_CDF()`, `plot_posterior_predictive()` and `get_probability_table()`
  provide posterior predictive distributions and probabilities for a single
  value or a range of values.
* `start_app_view()` launches the interactive Shiny application, which is called BSI (Bayesian Statistical Intervals), a different name from the package.
* The synthetic data set `data_InputFileTester2` is included for examples and
  tests.

## Changes

* The package is called `bayesintervals`; during development on GitHub it was called `BayesianStatisticalIntervalsCalculator` (BSI).
* The package is released under the Apache License 2.0.
* Dependencies were reduced: `ellipse`, `gridExtra`, `shinythemes` and
  `htmltools` are no longer used (confidence ellipses are computed
  internally and the application theme comes from `bslib`). `later` and
  `bslib` are now imported because the application needs them.
* Inputs are validated with clear error messages: interval levels, the multiplication factor, covariate arguments, column types, and the `by` column. Missing or infinite MCMC draws and negative variances now stop with an error; previously missing draws were silently treated as zero.
* Reducible imprecision is `NA`, with a warning, for a group where the tolerance bound equals the median (division by zero); it used to be `Inf` or `NaN` without a message. Duplicate `By` rows in the interval results now give a warning, and omitting `by` gives a message instead of a warning.* `calculate_predicted_distributions()` is much faster (about 15 times for several groups of 10,000 draws) and always returns numeric columns. Infinite covariate values are rejected.
* The documentation cites Krishnamoorthy and Mathew (2009) for the tolerance interval methods.
* Prediction intervals are now computed from the full posterior predictive distribution, the average of the per-draw normal distribution functions (manuscript Eq. 6-7), instead of one normal built from the averaged mean and variance. The new intervals are slightly wider when the draws disagree about the mean, i.e. in small samples; reducible imprecision changes accordingly. `calculate_PI()` uses a bracketed Newton-bisection search, and `Eta` is now a relative tolerance (default 1e-8) instead of a gap in probability.
* The native pipe `|>` replaces `%>%`.
* All examples are now runnable and use small synthetic data.
* `start_app_view()` restores the `shiny.maxRequestSize` option when the
  application is closed.
* The application no longer links to internal documentation.

## Bug fixes

* Two-sided Krishnamoorthy tolerance intervals (`use_krishnamoorthy = TRUE`):
  `calculate_gamma_difference()` now counts the draws whose own interval lies
  inside the candidate interval [g1, g2] (`lower_quantile >= g1` and
  `upper_quantile <= g2`), as in Krishnamoorthy and Mathew (2009). The
  comparison used to point the other way, which gave intervals that were too
  narrow.
* `find_optimal_point()` returned `gamma_achieved` as 100 minus the gap to the
  target, which is wrong for any target other than 100. It now reports the
  achieved percentage.
* `plot_krishnamoorthy()` and `quantile_scatterplot()` show the upper quantile
  on the x-axis and the lower quantile on the y-axis, so that the shaded
  region of draws inside the interval is the upper-left rectangle and the
  optimal point is at its lower-right corner. The orange line in
  `plot_krishnamoorthy()` now spans only the range of the draws, so the axes
  stay tight to the points at any scale of the data.
* The same changes were made in the Shiny application.
* `quantile_scatterplot()` now draws confidence ellipses when called from
  the package. The helper that computes them previously existed only inside
  the Shiny application, so the ellipses were silently skipped.
* Fixed functions that used `%>%` without importing it and would fail in a
  clean R session.
* Removed non-ASCII characters from R code and fixed broken
  cross-references in the documentation.
