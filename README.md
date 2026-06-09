# BayesianStatisticalIntervalsCalculator (BSI)

## Overview

BSI is a Bayesian statistical analysis tool designed for rapid calculation and visualization of statistical intervals from MCMC output data. It provides both an R package for programmatic access and a Shiny application for interactive analysis.

<img width="1833" height="779" alt="image" src="https://github.com/user-attachments/assets/94d58c42-2273-48ce-bf78-0587246ecce3" />


The package allows users to:
- Calculate and visualize credible, prediction, and tolerance intervals from MCMC samples
- Define Linear Mixed Models with fixed effects and variance components
- Support both normal and log-normal data transformations
- Generate interactive visualizations of statistical intervals
- Include error handling and data validation 

## Installation

The package is not yet available on CRAN. You can install it from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Sanofi-GitHub/BSI-APP/tree/master")
```

For access issues or questions about installation, please contact the development team or the maintainer:
Ekin Ören (Ekin.Oeren@sanofi.com)

## Quick Start

```r
library(BayesianStatisticalIntervalsCalculator)

# Example MCMC data
toy_data <- data.frame(
  EGTPTNUM = c(60, 60, 60, 60, 60, 90, 90, 90, 90, 90),
  s2 = c(74.18, 58.92, 75.04, 64.75, 67.81, 59.75, 51.91, 57.58, 64.76, 52.09),
  b2 = c(-5.23, -3.97, -2.92, -7.64, -2.57, 0.61, -2.02, -0.41, 0.08, 2.43),
  b4 = c(1.32, 3.68, 3.39, 3.13, 3.50, 5.57, 4.90, 1.71, 3.58, 3.36)
)

# Calculate intervals
result <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8
)

# View results
print(result$summary)
print(result$data)
```

## Main Features

### Statistical Interval Calculations

The core function `bayesian_statistical_intervals()` computes three types of intervals from MCMC samples:

- **Credible Intervals (CI)**: Bayesian confidence intervals for the median
- **Prediction Intervals (PI)**: Intervals for future individual observations
- **Tolerance Intervals (TI)**: Intervals containing a specified proportion of the population with given confidence

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | data.frame | - | MCMC samples dataset |
| `fixed_effects` | character vector | - | Column names for fixed effects |
| `random_params` | character vector | - | Column names for variance components |
| `by` | character | "Single Chain Dataset" | Grouping variable |
| `tolerance_interval_level` | numeric | - | TI probability level (0-100) |
| `credibility_of_tolerance_interval` | numeric | 95 | TI confidence level (0-100) |
| `credible_interval_level` | numeric | 95 | CI level (0-100) |
| `prediction_interval_level` | numeric | 90 | PI level (0-100) |
| `multiplication_factor` | numeric | 1 | Variance multiplier |
| `log_normal` | logical | FALSE | Apply log-normal transformation |
| `covariate_cols` | character vector | NULL | Covariate column names |
| `covariate_values` | numeric vector | NULL | Covariate multipliers |
| `use_krishnamoorthy` | logical | FALSE | Use Krishnamoorthy method for TI |
| `krishnamoorthy_tolerance` | numeric | 0.01 | Convergence tolerance |

**Example:**

```r
result <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  credible_interval_level = 90,
  prediction_interval_level = 90,
  multiplication_factor = 1.8
)

# Results include:
# - result$data: Data frame with all intervals
# - result$summary: Text summary of parameters
# - result$pred_dist: Predicted distributions
```

### Covariate Support

Covariates allow you to incorporate additional variables into your statistical model, adjusting predictions based on specific covariate values.

**Single Covariate:**

```r
result_with_cov <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8,
  covariate_cols = "b4",        # Add covariate
  covariate_values = 1.5        # Covariate multiplier
)
```

**Multiple Covariates:**

```r
result_multi_cov <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8,
  covariate_cols = c("b3", "b4"),      # Multiple covariates
  covariate_values = c(0.9, 1.5)      # Different multipliers
)
```

### Krishnamoorthy Two-Sided Tolerance Intervals

The Krishnamoorthy method provides an alternative approach to calculating two-sided tolerance intervals that may be more appropriate when the distribution of lower and upper quantiles is not symmetric.

```r
result_krishnamoorthy <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8,
  covariate_cols = "b4",
  covariate_values = 1.5,
  use_krishnamoorthy = TRUE  # Enable Krishnamoorthy method
)

# Results include K_TI_Lower and K_TI_Upper columns
print(result_krishnamoorthy$data)

# View Krishnamoorthy method details
print(result_krishnamoorthy$krishnamoorthy_results)
```

### Log-Normal Transformation

When working with data in log space, the package can automatically transform results back to the original scale.

```r
result_lognormal <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8,
  log_normal = TRUE  # Apply log-normal transformation
)
```

## Visualization Functions

### One-Sided Interval Plots

Create separate visualizations for upper and lower bounds of statistical intervals.

```r
# Calculate intervals
result <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  multiplication_factor = 1.8
)

# Create one-sided plots
one_sided_plots <- plot_one_sided_intervals(
  data = result$data,
  show_reducible = TRUE  # Show reducible uncertainty
)

# Display plots
one_sided_plots$upper
one_sided_plots$lower
```

### Quantile Scatterplot (Two-Sided Intervals)

Visualize the relationship between lower and upper quantiles with confidence ellipses.

```r
# Create ellipsoid scatterplot
ellipsoid_plot <- quantile_scatterplot(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_level = 90,
  ellipse_confidence = 80,
  multiplication_factor = 1.8,
  covariate_cols = "b4",
  covariate_values = 1.5
)

print(ellipsoid_plot)

# Or use pre-calculated distributions
ellipsoid_from_pred <- quantile_scatterplot(
  pred_dist = result$pred_dist,
  tolerance_level = 90,
  ellipse_confidence = 80
)
```

### Krishnamoorthy Plot

Visualize the Krishnamoorthy method results with scatter points, optimal line, and optimal point.

```r
# First calculate with Krishnamoorthy method
result_k <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  use_krishnamoorthy = TRUE
)

# Create Krishnamoorthy plot for specific group
krishna_plot <- plot_krishnamoorthy(
  data = result_k$pred_dist,
  by_value = "60",
  optimal_point = result_k$krishnamoorthy_results[["60"]],
  tolerance_level = 90,
  confidence_level = 80
)

print(krishna_plot)
```

### Posterior Predictive Distribution

Visualize the posterior predictive distribution with confidence intervals.

```r
# Create CDF data
cdf_data <- create_CDF(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_level = 90,
  n_points = 50,
  ci_level = 90,
  multiplication_factor = 1.8
)

# Basic plot
pp_plot <- plot_posterior_predictive(
  cdf_data = cdf_data[["60"]],
  title = "Posterior Predictive Distribution - EGTPTNUM 60"
)

# Plot with specific x-value
pp_plot_x <- plot_posterior_predictive(
  cdf_data = cdf_data[["60"]],
  x_value = "5",
  title = "Posterior Predictive Distribution with P(X ≤ 5)"
)

# Plot with probability range
pp_plot_range <- plot_posterior_predictive(
  cdf_data = cdf_data[["60"]],
  x_value = "[0, 10]",
  title = "Posterior Predictive Distribution with P(0 ≤ X ≤ 10)"
)
```

## Complete Workflow Example

Here's a complete analysis workflow demonstrating multiple features:

```r
library(BayesianStatisticalIntervalsCalculator)

# Step 1: Calculate all intervals with covariates and Krishnamoorthy method
complete_result <- bayesian_statistical_intervals(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_interval_level = 90,
  credibility_of_tolerance_interval = 80,
  credible_interval_level = 90,
  prediction_interval_level = 90,
  multiplication_factor = 1.8,
  covariate_cols = c("b3", "b4"),
  covariate_values = c(0.9, 1.5),
  use_krishnamoorthy = TRUE
)

# Step 2: View results
print(complete_result$summary)
print(complete_result$data)

# Step 3: Create one-sided interval plots
one_sided <- plot_one_sided_intervals(
  data = complete_result$data,
  show_reducible = TRUE
)
one_sided$upper
one_sided$lower

# Step 4: Create ellipsoid scatterplot
ellipsoid <- quantile_scatterplot(
  pred_dist = complete_result$pred_dist,
  tolerance_level = 90,
  ellipse_confidence = 80
)
print(ellipsoid)

# Step 5: Create Krishnamoorthy plots for each group
krishna_60 <- plot_krishnamoorthy(
  data = complete_result$pred_dist,
  by_value = "60",
  optimal_point = complete_result$krishnamoorthy_results[["60"]],
  tolerance_level = 90,
  confidence_level = 80
)
print(krishna_60)

# Step 6: Create posterior predictive distribution
cdf_data <- create_CDF(
  pred_dist = complete_result$pred_dist,
  n_points = 50,
  ci_level = 90
)

pp_plot <- plot_posterior_predictive(
  cdf_data = cdf_data[["60"]],
  x_value = "5"
)
print(pp_plot)
```

## Function Reference

| Function | Description |
|----------|-------------|
| `bayesian_statistical_intervals()` | Main function for calculating all statistical intervals |
| `calculate_predicted_distributions()` | Calculate distribution parameters from MCMC samples |
| `calculate_ci()` | Calculate credible intervals |
| `calculate_tolerance_interval()` | Calculate tolerance intervals |
| `calculate_PI()` | Calculate prediction intervals |
| `plot_one_sided_intervals()` | Create one-sided interval plots |
| `quantile_scatterplot()` | Create two-sided tolerance interval scatterplot |
| `plot_krishnamoorthy()` | Visualize Krishnamoorthy method results |
| `create_CDF()` | Generate CDF data for posterior predictive distribution |
| `plot_posterior_predictive()` | Plot posterior predictive distribution |
| `load_n_inspect_data()` | Load and validate data from files |
| `concatenate_by_columns()` | Merge multiple grouping columns |

## Data Requirements

The package expects MCMC output data in a data frame format with:

- **Numeric columns** for fixed effects and variance components
- **Grouping variable(s)** for stratified analysis (optional)
- **No missing values** in critical columns (warnings will be issued)

Example data structure:

```r
data.frame(
  Group = c("A", "A", "B", "B"),
  FixedEffect1 = c(1.2, 1.5, 2.1, 1.9),
  Variance1 = c(0.5, 0.6, 0.7, 0.8),
  Covariate1 = c(10, 12, 15, 14)
)
```

## Advanced Usage

### Pre-Calculated Distributions

For efficiency, you can calculate distributions once and reuse them:

```r
# Calculate distributions separately
pred_dist <- calculate_predicted_distributions(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  tolerance_level = 90,
  multiplication_factor = 1.8
)

# Use in multiple visualizations
ellipsoid <- quantile_scatterplot(pred_dist = pred_dist, tolerance_level = 90)
cdf_data <- create_CDF(pred_dist = pred_dist, n_points = 50)
```

### Filtering Specific Groups

Focus analysis on specific groups:

```r
# Show only specific groups in scatterplot
ellipsoid_filtered <- quantile_scatterplot(
  data = toy_data,
  fixed_effects = "b2",
  random_params = "s2",
  by = "EGTPTNUM",
  selected_by_values = c("60"),  # Show only group 60
  tolerance_level = 90,
  ellipse_confidence = 80
)
```

### Custom Visualization Options

```r
# Customize posterior predictive plot
pp_plot <- plot_posterior_predictive(
  cdf_data = cdf_data[["60"]],
  x_value = "[2, 8]",  # Range query
  show_ci = TRUE,
  line_colors = list(
    main = 'rgba(0,100,80,1)',
    ci = 'blue',
    ci_ribbon = 'rgba(0,0,255,0.2)'
  )
)
```

## Shiny Application

The package includes a companion Shiny application for interactive analysis:

```r
# Launch the Shiny app
start_app_view()
```

The Shiny app provides:
- Interactive data upload and preview
- Visual formula builder
- Real-time interval calculations
- Dynamic posterior predictive distribution plots
- Comprehensive error handling and validation

## Understanding MCMC Analysis and Statistical Intervals

### MCMC Samples

Markov Chain Monte Carlo (MCMC) methods generate samples from posterior distributions in Bayesian analysis. Each row in your dataset typically represents one sample from the posterior distribution of your model parameters.

### Statistical Intervals

- **Credible Intervals**: Summarize uncertainty about model parameters (analogous to frequentist confidence intervals)
- **Prediction Intervals**: Account for both parameter uncertainty and inherent variability in future observations
- **Tolerance Intervals**: Contain a specified proportion of the population with a given confidence level

### Linear Mixed Model

The package assumes a Linear Mixed Model structure:

```
y ~ N(fixed_effects + covariates, multiplication_factor * (random_effects))
```

Where:
- Fixed effects represent systematic effects
- Random effects represent variance components
- Covariates can be included with custom multipliers

## Citation

If you use this package in your research, please cite:

```
Ören, E. (2024). BayesianStatisticalIntervalsCalculator: 
Bayesian Statistical Intervals from MCMC Samples. 
R package version 0.x.x.
```

## License

This package is licensed under the [MIT License](LICENSE.md). Copyright © 2025 Sanofi. All rights reserved.

## Contact

For bug reports, feature requests, and other questions:
- GitHub Issues: https://github.com/Sanofi-Public/BayesianSI.Calculator/issues
- Maintainer: Ekin Ören (Ekin.Oeren@sanofi.com)

---

**Coming Soon:**

- [Internal_BSI Shiny App]
- [Technical Documentation]
