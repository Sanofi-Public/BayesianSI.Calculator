#' Create Krishnamoorthy Two-Sided Tolerance Interval Plot
#'
#' Creates an interactive plot showing the Krishnamoorthy method results for two-sided
#' tolerance intervals, including scatter points, optimal line, and optimal point.
#'
#' @param data A data frame containing lower and upper quantiles
#' @param by_value Character or numeric, the specific by-value to plot
#' @param optimal_point List containing Krishnamoorthy method results with elements:
#'        \itemize{
#'          \item g1: optimal lower bound
#'          \item g2: optimal upper bound
#'          \item gamma_achieved: achieved confidence level
#'        }
#' @param tolerance_level Numeric, tolerance level percentage (between 0 and 100)
#' @param confidence_level Numeric, confidence level percentage (between 0 and 100)
#' @param log_normal Logical, whether to transform results from log space (default = FALSE)
#'
#' @return A plotly object containing the Krishnamoorthy plot
#'
#' @references
#' Krishnamoorthy K, Mathew T (2009). \emph{Statistical Tolerance Regions: Theory,
#' Applications, and Computation}. Wiley, Hoboken, NJ. \doi{10.1002/9780470473900}.
#'
#' @examples
#' set.seed(1)
#' mcmc <- data.frame(
#'   Group = rep(c("A", "B"), each = 50),
#'   beta = c(rnorm(50, 10, 0.5), rnorm(50, 12, 0.5)),
#'   sigma2 = rgamma(100, shape = 20, rate = 10),
#'   age = rnorm(100, 50, 5)
#' )
#'
#' pred_dist <- calculate_predicted_distributions(
#'   mcmc, fixed_effects = "beta", random_params = "sigma2",
#'   by = "Group", tolerance_level = 90
#' )
#' gP_mean <- calculate_gP_mean(pred_dist, 90)
#' optimal <- find_optimal_point(pred_dist, gP_mean, gamma_actual = 80)
#' p <- plot_krishnamoorthy(
#'   data = pred_dist, by_value = "A", optimal_point = optimal[["A"]],
#'   tolerance_level = 90, confidence_level = 80
#' )
#'
#' @importFrom plotly plot_ly add_trace layout
#' @export
plot_krishnamoorthy <- function(data, by_value, optimal_point,
                                tolerance_level, confidence_level,
                                log_normal = FALSE) {

  # Input validation
  if (!all(c("lower_quantile", "upper_quantile", "By") %in% names(data))) {
    stop("Data must contain 'lower_quantile', 'upper_quantile', and 'By' columns")
  }

  if (is.null(optimal_point) || !all(c("g1", "g2", "gamma_achieved") %in% names(optimal_point))) {
    stop("optimal_point must contain 'g1', 'g2', and 'gamma_achieved'")
  }

  # Filter data for the selected by_value
  data_subset <- data[data$By == by_value, ]

  if (nrow(data_subset) == 0) {
    stop(paste("No data found for by_value:", by_value))
  }

  # Calculate percentiles for labels
  lower_percentile <- (100 - tolerance_level) / 2
  upper_percentile <- 100 - lower_percentile

  # Initialize plot
  p <- plot_ly()

  # Get optimal point values
  g1_log <- if(log_normal) log(optimal_point$g1) else optimal_point$g1
  g2_log <- if(log_normal) log(optimal_point$g2) else optimal_point$g2
  gP_mean_log <- (g1_log + g2_log) / 2

  # Determine display values
  g1_display <- if(log_normal) exp(g1_log) else g1_log
  g2_display <- if(log_normal) exp(g2_log) else g2_log

  # Add scatter points
  p <- add_trace(p,
                 data = data_subset,
                 x = if(log_normal) ~exp(upper_quantile) else ~upper_quantile,
                 y = if(log_normal) ~exp(lower_quantile) else ~lower_quantile,
                 type = 'scatter',
                 mode = 'markers',
                 name = paste("By:", by_value),
                 marker = list(
                   size = 8,
                   color = '#1f77b4',
                   opacity = 0.8
                 ),
                 hoverinfo = 'text',
                 text = ~paste(
                   "By:", by_value,
                   "<br>Lower Quantile:", round(if(log_normal) exp(lower_quantile) else lower_quantile, 4),
                   "<br>Upper Quantile:", round(if(log_normal) exp(upper_quantile) else upper_quantile, 4)
                 )
  )

  # Calculate line points
  x_range <- range(data_subset$upper_quantile)
  y_range <- range(data_subset$lower_quantile)
  # The line spans the padded data range only, so the axes stay tight to the point cloud
  # whatever the scale of the data
  x_line_pad <- (x_range[2] - x_range[1]) * 0.1
  y_line_pad <- (y_range[2] - y_range[1]) * 0.1
  x_line_log <- seq(x_range[1] - x_line_pad, x_range[2] + x_line_pad, length.out = 400)
  y_line_log <- 2 * gP_mean_log - x_line_log
  in_box <- y_line_log >= y_range[1] - y_line_pad & y_line_log <= y_range[2] + y_line_pad
  if (any(in_box)) x_line_log <- x_line_log[in_box]

  if(log_normal) {
    x_line_display <- exp(x_line_log)
    y_line_display <- exp(-log(x_line_display) + 2*gP_mean_log)
  } else {
    x_line_display <- x_line_log
    y_line_display <- 2 * gP_mean_log - x_line_log
  }

  # Add Krishnamoorthy line
  p <- add_trace(p,
                 x = x_line_display,
                 y = y_line_display,
                 type = 'scatter',
                 mode = 'lines',
                 name = 'Krishnamoorthy Line',
                 line = list(
                   color = '#FF8C00',
                   width = 2,
                   dash = 'solid'
                 ),
                 hoverinfo = 'text',
                 text = "Krishnamoorthy Line"
  )

  # Add optimal point
  p <- add_trace(p,
                 x = g2_display,
                 y = g1_display,
                 type = 'scatter',
                 mode = 'markers',
                 name = 'Optimal Point',
                 marker = list(
                   size = 10,
                   color = '#FF8C00',
                   symbol = 'diamond'
                 ),
                 hoverinfo = 'text',
                 text = sprintf(
                   "Two-sided TI (Krishnamoorthy)<br>Lower TI Limit: %.4f<br>Upper TI Limit: %.4f<br>Achieved Credible Level: %.2f%%",
                   g1_display, g2_display, optimal_point$gamma_achieved
                 )
  )

  # Calculate display ranges and padding
  x_min_display <- if(log_normal) exp(min(x_range)) else min(x_range)
  x_max_display <- if(log_normal) exp(max(x_range)) else max(x_range)
  y_min_display <- if(log_normal) exp(min(y_range)) else min(y_range)
  y_max_display <- if(log_normal) exp(max(y_range)) else max(y_range)

  x_padding <- (x_max_display - x_min_display) * 0.1
  y_padding <- (y_max_display - y_min_display) * 0.1

  # Set layout with shapes
  p <- layout(p,
              title = sprintf("Two-Sided Tolerance Intervals (%.1f%% Probability Level, %.1f%% Credible Level)",
                              tolerance_level, confidence_level),
              xaxis = list(
                title = sprintf("Upper Tolerance Interval Probability Level (%.1f%%)%s",
                                upper_percentile, if(log_normal) " [Exponentiated]" else ""),
                zeroline = TRUE
              ),
              yaxis = list(
                title = sprintf("Lower Tolerance Interval Probability Level (%.1f%%)%s",
                                lower_percentile, if(log_normal) " [Exponentiated]" else ""),
                zeroline = TRUE
              ),
              hovermode = 'closest',
              showlegend = TRUE,
              legend = list(
                orientation = "h",
                xanchor = "center",
                x = 0.5,
                y = -0.2
              ),
              shapes = list(
                # Vertical line at g2_display
                list(
                  type = "line",
                  x0 = g2_display,
                  x1 = g2_display,
                  y0 = y_min_display - y_padding,
                  y1 = y_max_display + y_padding,
                  line = list(
                    color = '#006400',
                    width = 2,
                    dash = 'dash'
                  )
                ),
                # Horizontal line at g1_display
                list(
                  type = "line",
                  x0 = x_min_display - x_padding,
                  x1 = x_max_display + x_padding,
                  y0 = g1_display,
                  y1 = g1_display,
                  line = list(
                    color = '#006400',
                    width = 2,
                    dash = 'dash'
                  )
                ),
                # Shaded area (upper left)
                list(
                  type = "rect",
                  x0 = x_min_display - x_padding,
                  x1 = g2_display,
                  y0 = g1_display,
                  y1 = y_max_display + y_padding,
                  fillcolor = "rgba(0, 100, 0, 0.1)",
                  line = list(width = 0)
                )
              )
  )

  return(p)
}
