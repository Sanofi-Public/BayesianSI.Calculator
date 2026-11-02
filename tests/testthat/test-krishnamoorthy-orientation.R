# Orientation of the two-sided tolerance interval construction (Krishnamoorthy and
# Mathew, 2009): a draw counts towards the credibility if its own interval
# [lower_quantile, upper_quantile] lies inside the candidate interval [g1, g2].
# In the plots the upper quantile is on the horizontal axis and the lower quantile on
# the vertical axis, so the shaded region (upper quantile <= g2, lower quantile >= g1)
# is the upper-left rectangle.

test_that("calculate_gamma_difference counts draws whose interval lies inside [g1, g2]", {
  sd <- data.frame(By = "A",
                   lower_quantile = c(-1, -1, -1, -3),
                   upper_quantile = c(1, 1, 1, 3))
  # g1 = -1.5, g2 = 1.5: the first three draws are inside, the fourth is not
  res <- calculate_gamma_difference(sd, gP_mean = 0, point_x_vector = -1.5, gamma_actual = 75)
  expect_equal(res$points_satisfying, 3)
  expect_equal(unname(res$gamma_differences), 0)
  # a very wide candidate contains every draw, a very narrow one none
  expect_equal(calculate_gamma_difference(sd, 0, -10, 75)$points_satisfying, 4)
  expect_equal(calculate_gamma_difference(sd, 0, -0.5, 75)$points_satisfying, 0)
})

test_that("the share of draws inside the optimal interval matches the target credibility", {
  set.seed(42)
  mcmc <- data.frame(g = "A", b = rnorm(400, 5, 0.6), s2 = rgamma(400, 30, 15))
  pd <- calculate_predicted_distributions(mcmc, "b", "s2", "g", tolerance_level = 90)
  gp <- calculate_gP_mean(pd, 90)
  opt <- find_optimal_point(pd, gp, gamma_actual = 80)[["A"]]
  inside <- mean(pd$lower_quantile >= opt$g1 & pd$upper_quantile <= opt$g2)
  expect_equal(inside, 0.80, tolerance = 0.01)
  expect_equal(unname((opt$g1 + opt$g2) / 2), unname(gp["A"]), tolerance = 1e-8)
  expect_lt(opt$g1, opt$g2)
})

test_that("gamma_achieved reports the share of draws inside the interval, not 100 minus the gap", {
  set.seed(42)
  mcmc <- data.frame(g = "A", b = rnorm(400, 5, 0.6), s2 = rgamma(400, 30, 15))
  pd <- calculate_predicted_distributions(mcmc, "b", "s2", "g", tolerance_level = 90)
  gp <- calculate_gP_mean(pd, 90)
  opt <- find_optimal_point(pd, gp, gamma_actual = 80)[["A"]]
  expect_equal(unname(opt$gamma_achieved), unname(100 * opt$points_satisfying / opt$total_points))
  expect_equal(unname(opt$gamma_achieved), 80, tolerance = 0.01)
})

test_that("the Krishnamoorthy line stays within the padded data range at any data scale", {
  set.seed(11)
  for (centre in c(0, 5, 1500)) {
    mcmc <- data.frame(g = "A", b = rnorm(300, centre, 20), s2 = rgamma(300, 20, 0.01))
    pd <- calculate_predicted_distributions(mcmc, "b", "s2", "g", tolerance_level = 90)
    gp <- calculate_gP_mean(pd, 90)
    opt <- find_optimal_point(pd, gp, gamma_actual = 80)[["A"]]
    pk <- plotly::plotly_build(plot_krishnamoorthy(pd, "A", opt, 90, 80))
    ln <- Filter(function(tr) identical(tr$name, "Krishnamoorthy Line"), pk$x$data)[[1]]
    xr <- range(pd$upper_quantile); yr <- range(pd$lower_quantile)
    lx <- as.numeric(ln$x); ly <- as.numeric(ln$y)
    expect_gte(min(lx), xr[1] - 0.1 * diff(xr) - 1e-8)
    expect_lte(max(lx), xr[2] + 0.1 * diff(xr) + 1e-8)
    expect_gte(min(ly), yr[1] - 0.1 * diff(yr) - 1e-8)
    expect_lte(max(ly), yr[2] + 0.1 * diff(yr) + 1e-8)
    expect_equal(lx + ly, rep(unname(opt$g1 + opt$g2), length(lx)), tolerance = 1e-8)  # still the line g1 + g2 = const
  }
})

axis_title <- function(ax) if (is.list(ax$title)) ax$title$text else ax$title

test_that("plots have the upper quantile on the x-axis and the lower quantile on the y-axis", {
  set.seed(7)
  mcmc <- data.frame(g = rep(c("A", "B"), each = 60),
                     b = c(rnorm(60, 5, 0.5), rnorm(60, 7, 0.5)), s2 = rgamma(120, 30, 15))
  pd <- calculate_predicted_distributions(mcmc, "b", "s2", "g", tolerance_level = 90)
  gp <- calculate_gP_mean(pd, 90)
  opt <- find_optimal_point(pd, gp, gamma_actual = 80)[["A"]]
  A <- pd[pd$By == "A", ]

  # Krishnamoorthy plot
  pk <- plotly::plotly_build(plot_krishnamoorthy(pd, "A", opt, 90, 80))
  scatter <- pk$x$data[[1]]
  expect_equal(as.numeric(scatter$x), A$upper_quantile)
  expect_equal(as.numeric(scatter$y), A$lower_quantile)
  diamond <- Filter(function(tr) identical(tr$name, "Optimal Point"), pk$x$data)[[1]]
  expect_equal(as.numeric(diamond$x), unname(opt$g2))
  expect_equal(as.numeric(diamond$y), unname(opt$g1))
  expect_match(axis_title(pk$x$layout$xaxis), "^Upper Tolerance Interval")
  expect_match(axis_title(pk$x$layout$yaxis), "^Lower Tolerance Interval")
  rect <- Filter(function(s) identical(s$type, "rect"), pk$x$layout$shapes)[[1]]
  expect_equal(as.numeric(rect$x1), as.numeric(opt$g2))   # right edge of the shaded region at g2
  expect_equal(as.numeric(rect$y0), as.numeric(opt$g1))   # lower edge at g1; it extends to the top (upper left)
  expect_gt(rect$y1, rect$y0)
  expect_lt(rect$x0, rect$x1)

  # The shaded region holds the draws that the credibility refers to
  in_region <- A$upper_quantile <= rect$x1 & A$lower_quantile >= rect$y0
  expect_equal(mean(in_region), unname(opt$points_satisfying) / opt$total_points)

  # Quantile scatterplot with ellipses
  ps <- plotly::plotly_build(quantile_scatterplot(pred_dist = pd, tolerance_level = 90,
                                                  ellipse_confidence = 90))
  expect_equal(as.numeric(ps$x$data[[1]]$x), A$upper_quantile)
  expect_equal(as.numeric(ps$x$data[[1]]$y), A$lower_quantile)
  expect_match(axis_title(ps$x$layout$xaxis), "^Upper Tolerance Interval")
  expect_match(axis_title(ps$x$layout$yaxis), "^Lower Tolerance Interval")
  mean_tr <- Filter(function(tr) identical(tr$name, "Mean"), ps$x$data)[[1]]
  expect_equal(as.numeric(mean_tr$x), mean(A$upper_quantile))
  expect_equal(as.numeric(mean_tr$y), mean(A$lower_quantile))
})
