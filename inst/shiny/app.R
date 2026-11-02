# Bayesian Statistical Intervals (BSI) Shiny application
# This Shiny app allows users to upload CSV data, perform statistical calculations,
# and visualize the results.

# Index:
# 1. Library Imports
# 2. Global Variables and Functions
# 3. UI Definition
#    3.1 Main Layout
#    3.2 Sidebar Panel
#    3.3 Main Panel
# 4. Server Logic
#    4.1 Data Loading and Inspection
#    4.2 Reactive Values
#    4.3 Observe events
#    4.4 Outputs
# 5. Shiny App Execution


### SETUP ###
library(shiny)
library(bayesintervals)
library(DT)
library(bslib)
library(dplyr)

#for reading in large files
library(readr)
library(data.table)

#for plotting needs
library(ggplot2)
library(plotly)

# Computes confidence ellipse coordinates from covariance and center without external dependencies.
create_confidence_ellipse_points <- function(cov_matrix,
                                             centre,
                                             level = 0.95,
                                             npoints = 100) {
  if (!is.matrix(cov_matrix) || any(dim(cov_matrix) != c(2, 2))) {
    stop("'cov_matrix' must be a 2x2 matrix")
  }
  if (!is.numeric(centre) || length(centre) != 2) {
    stop("'centre' must be a numeric vector of length 2")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("'level' must be a number between 0 and 1")
  }
  if (!is.numeric(npoints) || length(npoints) != 1 || npoints < 4) {
    stop("'npoints' must be a number >= 4")
  }

  eig <- eigen(cov_matrix, symmetric = TRUE)
  if (any(!is.finite(eig$values)) || any(eig$values <= .Machine$double.eps)) {
    stop("Covariance matrix must be positive definite")
  }

  radius <- sqrt(stats::qchisq(level, df = 2))
  angles <- seq(0, 2 * pi, length.out = npoints)
  unit_circle <- rbind(cos(angles), sin(angles))

  scale_matrix <- eig$vectors %*% diag(sqrt(eig$values), nrow = 2)
  transformed <- radius * (scale_matrix %*% unit_circle)

  points <- cbind(centre[1] + transformed[1, ], centre[2] + transformed[2, ])
  colnames(points) <- c("x", "y")
  points
}


# Input checks. The package keeps an identical copy in R/validate_inputs.R so
# that this application can run on its own; keep the two in sync.
check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a single positive number", name), call. = FALSE)
  }
  invisible(TRUE)
}

check_draw_columns <- function(dt, cols, label) {
  if (nrow(dt) == 0) {
    stop("'data' has no rows", call. = FALSE)
  }
  for (col in cols) {
    x <- dt[[col]]
    if (!is.numeric(x)) {
      converted <- suppressWarnings(as.numeric(as.character(x)))
      if (any(is.na(converted) & !is.na(x))) {
        stop(sprintf("%s column '%s' must be numeric but contains non-numeric values", label, col),
             call. = FALSE)
      }
      x <- converted
    }
    bad <- which(!is.finite(x))
    if (length(bad) > 0) {
      stop(sprintf("%s column '%s' contains %d missing or non-finite value(s) (first at row %d)",
                   label, col, length(bad), bad[1]), call. = FALSE)
    }
  }
  invisible(TRUE)
}

check_variance <- function(variance_vector) {
  neg <- which(variance_vector < 0)
  if (length(neg) > 0) {
    stop(sprintf("The variance components sum to a negative value in %d draw(s) (first at row %d); variances must be non-negative",
                 length(neg), neg[1]), call. = FALSE)
  }
  invisible(TRUE)
}

# Reducible imprecision helpers (copy of R/validate_inputs.R; keep in sync)
reducible_ratio <- function(bound, pi, median_value, side, warn = TRUE) {
  ratio <- (bound - pi) / (bound - median_value)
  undefined <- is.finite(bound) & is.finite(pi) & is.finite(median_value) & !is.finite(ratio)
  if (any(undefined)) {
    ratio[undefined] <- NA_real_
    if (warn) {
      warning(sprintf(
        "Reducible imprecision (%s) is not defined for %d group(s) because the tolerance bound equals the median; NA is returned",
        side, sum(undefined)), call. = FALSE)
    }
  }
  ratio
}

add_reducible_imprecision <- function(result, warn = TRUE, percent = FALSE) {
  upper <- reducible_ratio(result$TI_Upper, result$PI_Upper, result$Median, "upper", warn)
  lower <- reducible_ratio(result$TI_Lower, result$PI_Lower, result$Median, "lower", warn)
  if (percent) {
    result$ReducibleUpper <- round(upper * 100, 3)
    result$ReducibleLower <- round(lower * 100, 3)
  } else {
    result$ReducibleUpper <- round(upper, 5)
    result$ReducibleLower <- round(lower, 5)
  }
  result
}

### UI ###
ui <- fluidPage(
  theme = bs_theme(preset = "darkly"),

  tags$head(

    tags$style(HTML("
      .btn-primary {
        background-color: #668167;
        border-color: #668167;
        margin-top: 10px;
        margin-bottom: 10px;
      }
      .btn-primary:hover, .btn-primary:focus, .btn-primary:active {
        background-color: #7a9879;
        border-color: #7a9879;
      }
      .js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar,
      .js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar,
      .js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar,
      .js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {
        background: #4a6b4a;
        border-color: #4a6b4a;
      }
      #formula_output {
        color: white;
        font-weight: bold;
        font-size: 18px;
      }
      #formula_output .function {
        color: #90EE90;
      }
      .dataTables_wrapper .dataTables_info {
        color: white;
      }
      .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter {
        color: white;
      }
      .dataTables_wrapper .dataTables_length select, .dataTables_wrapper .dataTables_filter input {
        color: black;
      }
      table.dataTable thead th {
        text-align: center !important;
      }
      table.dataTable tbody td {
        text-align: center !important;
      }
      table.dataTable tbody td:first-child {
        text-align: left !important;
      }
      .load-data-container {
        display: flex;
        align-items: center;
        justify-content: space-between;
        margin-bottom: 15px;
      }
      .log-normal-container {
        display: flex;
        align-items: center;
        justify-content: flex-start;
      }
      .checkbox {
        margin: 0;
        padding: 0;
      }
      .checkbox label {
        display: flex;
        align-items: center;
        margin: 0;
        padding: 0;
      }
      #log_normal, #show_reducible {
        -webkit-appearance: none;
        -moz-appearance: none;
        appearance: none;
        width: 14px;
        height: 14px;
        border: 2px solid #668167;
        border-radius: 4px;
        outline: none;
        transition: background-color 0.3s;
        vertical-align: middle;
        margin: 0;
        padding: 0;
      }
      #log_normal:checked, #show_reducible:checked {
        background-color: #668167;
      }
      #log_normal:checked::after, #show_reducible:checked::after {
        content: '';
        display: block;
        width: 3px;
        height: 7px;
        border: solid white;
        border-width: 0 2px 2px 0;
        transform: rotate(45deg);
        margin: 1px 0 0 4px;
      }
      #log_normal + span, #show_reducible + span {
        margin-left: 5px;
        line-height: 14px;
      }

      /* New CSS rules added for smaller text and compact layout in distribution plots controls */
      .small-text-controls .form-group * {
        font-size: 0.9em;
      }
      .small-text-controls .control-label {
        font-size: 0.9em;
        font-weight: normal;
        margin-bottom: 2px;
      }
      .small-text-controls .form-group {
        margin-bottom: 5px;
      }

      /* New CSS rules added for margin under the intervals table */
      #ci_output {
        margin-bottom: 30px;
      }

      /* New CSS rules for quantile plot controls */

      .well {
        border-radius: 4px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.12);
      }

      /* Checkbox styling */
      .checkbox {
        margin: 0;
        padding: 0;
      }

      .checkbox label {
        display: flex;
        align-items: center;
        color: #ffffff;
        font-weight: normal;
      }

      .checkbox input[type='checkbox'] {
        margin-right: 8px;
      }

      /* Quantile scatterplot controls styling */
      .checkbox-container {
        display: flex;
        flex-wrap: wrap;
        gap: 15px;
        margin-bottom: 25px;
        padding: 10px;
        background-color: #000000;
        border-radius: 4px;
      }

      .checkbox-item {
        background-color: #000000;
        padding: 8px 15px;
        border-radius: 4px;
        border: 1px solid #668167;
      }

      .slider-section {
        margin-top: 25px;
        padding: 15px;
        background-color: #000000;
        border-radius: 4px;
      }

      .slider-description {
        color: #ffffff;
        margin-top: 10px;
        font-size: 0.9em;
      }

      .info-section {
    margin-top: 25px;
    padding: 15px;
    background-color: #000000;
    border-radius: 4px;
    border-left: 4px solid #668167;
    }

    .info-section p {
        margin: 0;
        font-size: 0.9em;
        line-height: 1.5;
    }

    "))
  ),


  div(
    style = "margin-bottom: 15px;",  # Add some space below the title section
    h1("Bayesian Statistical Intervals (BSI)", style = "margin-bottom: 5px;")  # Main title with reduced bottom margin
  ),

  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      div(
        class = "load-data-container",
        div(
          style = "flex: 1;",
          actionButton("load_data", "Load Data", class = "btn-primary")
        ),
        div(
          class = "log-normal-container",
          style = "flex: 1; display: flex; justify-content: flex-start; align-items: center; margin-left: 20px;",
          div(
            style = "display: flex; align-items: center;",
            tags$input(type = "checkbox", id = "log_normal"),
            tags$label(
              "for" = "log_normal",
              style = "margin-left: 5px; margin-bottom: 0;",
              "Log Normal"
            )
          )
        )
      ),

      selectInput("fixed_var", "Fixed Effects", choices = NULL, multiple = TRUE),

      fluidRow(
        column(6,
               selectInput("covariate_1", "Covariate",
                           choices = NULL,
                           multiple = FALSE)
        ),
        column(6,
               numericInput("covariate_value_1", "Covariate Value",
                            value = NA,
                            step = 0.1)
        )
      ),
      # Covariate Pair 2
      fluidRow(
        column(6,
               selectInput("covariate_2", "Covariate",
                           choices = NULL,
                           multiple = FALSE)
        ),
        column(6,
               numericInput("covariate_value_2", "Covariate Value",
                            value = NA,
                            step = 0.1)
        )
      ),
      # Covariate Pair 3
      fluidRow(
        column(6,
               selectInput("covariate_3", "Covariate",
                           choices = NULL,
                           multiple = FALSE)
        ),
        column(6,
               numericInput("covariate_value_3", "Covariate Value",
                            value = NA,
                            step = 0.1)
        )
      ),

      selectInput("random_var", "Variances", choices = NULL, multiple = TRUE),

      numericInput("multiplication_factor", "Multiplication Factor for Variances",
                   value = 1, min = 0, step = 0.1),

      selectInput("by", "By Identifier", choices = c("Single Chain Dataset"), selected = "Single Chain Dataset", multiple = TRUE),
      uiOutput("byHelpText"),

      sliderInput("tolerance_level",
                  HTML("Tolerance Interval (%)<br>Probability Level"),
                  min = 0, max = 100, value = 90, step = 0.5),

      sliderInput("confidence_of_tolerance",
                  HTML("Tolerance Interval (%)<br>Credibility Level"),
                  min = 0, max = 100, value = 80, step = 0.5),

      sliderInput("percent_for_pi", "Prediction Interval (%)",
                  min = 0, max = 100, value = 90, step = 0.5),

      sliderInput("percent_for_ci", "Credible Interval (%)",
                  min = 0, max = 100, value = 90, step = 0.5),

      actionButton("calculate_ci", "Calculate Intervals", class = "btn-primary"),

      div(
        class = "log-normal-container",
        style = "flex: 1; display: flex; justify-content: flex-start; align-items: center; margin-left: 10px;",
        div(
          style = "display: flex; align-items:center;",
          tags$input(type="checkbox", id = "show_reducible"),
          tags$label(
            "for" = "show_reducible",
            style = "margin-left: 5px; margin-bottom: 0px;",
            "Reducible Imprecision"
          )
        )
      ),

      div(
        style = "display: flex; flex-direction: column;",
        actionButton("plot_quantiles", "Two-Sided Tolerance Intervals", class = "btn-primary"),
        actionButton("plot_distribution", "Posterior Predictive Distribution", class = "btn-primary")
      )


    ),

    mainPanel(

      tabsetPanel(id = "main_tabs",
                  tabPanel("Data Preview", DTOutput("data_preview")),
                  tabPanel("Formula", uiOutput("formula_output")),
                  tabPanel("One-Sided Intervals",
                           DTOutput("ci_output"),
                           plotlyOutput("intervals_plot_upper"),
                           div(style = "margin: 20px 0;"),  # Adds 20px margin above and below to seperate the two internal plots
                           plotlyOutput("intervals_plot_lower"),
                           tags$div(
                             style = "margin-top: 20px; padding: 15px; background-color: #344b61; border-radius: 4px; color: #ffffff;",
                             tags$p(
                               HTML("For two-sided tolerance intervals please go to the <b>Two-Sided Tolerance Intervals</b> tab.")
                             )
                           )
                  ),
                  tabPanel("Two-Sided Intervals",
                           fluidRow(
                             column(12,
                                    h4("Tolerance Interval Scatterplot Ellipsoid Controls"),
                                    uiOutput("quantile_scatter_controls")
                             )
                           ),
                           plotlyOutput("quantile_scatter_plot", height = "600px"),
                           fluidRow(
                             column(12,
                                    h4("Two-Sided Tolerance Interval Krishnamoorthy Method Scatterplot Controls"),
                                    uiOutput("krishnamoorthy_controls"),
                             )
                           ),
                           plotlyOutput("krishnamoorthy_plot", height = "600px"),
                           fluidRow(
                             column(12,
                                    h4("Two-Sided Intervals"),
                                    DTOutput("quantile_ci_output")
                             )
                           )
                  ),
                  tabPanel("Post. Pred. Distribution",
                           uiOutput("distribution_plots")),
                  tabPanel("Warning & Errors", verbatimTextOutput("error_output"))
      ),

      # Add the results table output outside the tabset panel but with display:none
      div(style = "display: none;", DTOutput("results_table"))
    )
  )
)



### SERVER ###
server <- function(input, output, session) {


  ### Variables
  InputData <- reactiveVal()
  ErrorMessage <- reactiveVal("")
  IntervalResults <- reactiveVal()
  FormulaGenerated <- reactiveVal(FALSE)
  CICalculated <- reactiveVal(FALSE)
  firstInputGiven <- reactiveVal(FALSE)
  ConfidenceOfTolerance <- reactiveVal(95)
  MultiplicationFactorVariances <- reactiveVal(1)
  DataInspectionResults <- reactiveVal(NULL)
  predictedDistributions <- reactiveVal(NULL)
  plot_data <- reactiveVal(NULL)
  intervalCalculationInitiated <- reactiveVal(FALSE)
  IntervalsOutputInitial <- reactiveVal(NULL)
  IntervalsOutputActual <- reactiveVal(NULL)
  byIdentifierChanged <- reactiveVal(FALSE)
  # Reactive time variables
  input_changed <- reactiveVal(FALSE)
  last_change_time <- reactiveVal(Sys.time())
  update_timer <- reactiveTimer(1000)
  # secondary reactives
  formulaString <- reactive({
    req(input$fixed_var, input$random_var, input$multiplication_factor)
    req(!reactivity_paused())

    print("Updating formula string")

    tryCatch({
      # Build covariate terms
      covariate_terms <- c()
      for(i in 1:3) {
        covariate_col <- input[[paste0("covariate_", i)]]
        covariate_val <- input[[paste0("covariate_value_", i)]]

        if(!is.null(covariate_col) && covariate_col != "None") {
          multiplier <- if(is.null(covariate_val) || is.na(covariate_val)) 1 else covariate_val
          if(multiplier == 1) {
            covariate_terms <- c(covariate_terms, covariate_col)
          } else {
            covariate_terms <- c(covariate_terms, paste0(multiplier, " * ", covariate_col))
          }
        }
      }

      covariate_string <- if(length(covariate_terms) > 0) {
        paste(covariate_terms, collapse = " + ")
      } else {
        ""
      }

      fixed_effects <- paste(input$fixed_var, collapse = " + ")
      random_effects <- paste(input$random_var, collapse = " + ")

      mult_factor <- ifelse(input$multiplication_factor != 1,
                            paste0(input$multiplication_factor, " * "),
                            "")

      # Build the mean part of the formula
      mean_parts <- c()
      if(covariate_string != "") mean_parts <- c(mean_parts, covariate_string)
      if(fixed_effects != "") mean_parts <- c(mean_parts, fixed_effects)

      mean_part <- if(length(mean_parts) > 0) {
        paste(mean_parts, collapse = " + ")
      } else {
        "<span style='color:red;'>Fixed Effects</span>"
      }

      formula <- paste0(
        "Calculating Probability Distributions by:<br><br>",
        "y ~ N(",
        mean_part,
        " , ",
        mult_factor,
        "(",
        ifelse(random_effects == "",
               "<span style='color:red;'>Random Effects</span>",
               random_effects),
        "))"
      )

      print("Generated formula:")
      print(formula)

      return(HTML(formula))
    }, error = function(e) {
      print(paste("Error in formulaString:", e$message))
      return(HTML("Error generating formula. Please check your inputs."))
    })
  })

  ForCDF1Dataframe <- reactive({
    req(predictedDistributions(), input$plot_distribution)
    req(!reactivity_paused())

    by_values <- unique(predictedDistributions()$By)

    result <- lapply(by_values, function(by_value) {
      if (isTRUE(input[[paste0("show_", by_value)]])) {
        data <- predictedDistributions()[predictedDistributions()$By == by_value, ]
        n_points <- input[[paste0("points_", by_value)]]

        # Validate number of points
        if (is.na(n_points) || !is.numeric(n_points) || n_points < 3 || n_points > 300) {
          return(NULL)
        }

        ci <- input[[paste0("ci_", by_value)]] / 100

        # Use mean of means for centering
        overall_mean <- mean(data$mean)
        overall_sd <- sqrt(mean(data$variance))

        # Generate x values
        xmin_log <- qnorm(0.0001, mean = overall_mean, sd = overall_sd)
        xmax_log <- qnorm(0.9999, mean = overall_mean, sd = overall_sd)
        x_values_log <- seq(xmin_log, xmax_log, length.out = n_points)

        # Transform x values for display if log_normal is TRUE
        x_values_display <- if(input$log_normal) exp(x_values_log) else x_values_log
        xmin_display <- if(input$log_normal) exp(xmin_log) else xmin_log
        xmax_display <- if(input$log_normal) exp(xmax_log) else xmax_log

        # Calculate CDF for each posterior sample at each x-value
        y_values <- do.call(rbind, lapply(1:nrow(data), function(i) {
          pnorm(x_values_log, mean = data$mean[i], sd = sqrt(data$variance[i]))
        }))

        # KEY CHANGE: Use median instead of mean for central tendency
        y_median <- apply(y_values, 2, median)

        # Calculate credibility intervals
        if (nrow(y_values) > 1) {
          alpha <- 1 - ci
          # Use standard quantile method (with interpolation)
          empirical_quantiles <- apply(y_values, 2, function(col) {
            quantile(col, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
          })
          ci_lower <- empirical_quantiles[1,]
          ci_upper <- empirical_quantiles[2,]
        } else {
          # Only one sample - no uncertainty
          ci_lower <- y_median
          ci_upper <- y_median
        }

        # Create result data frame
        result <- data.frame(
          By = by_value,
          x = x_values_display,
          x_log = x_values_log,
          y_mean = y_median,  # Note: Keep column name for compatibility
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          xmin = xmin_display,
          xmax = xmax_display,
          xmin_log = xmin_log,
          xmax_log = xmax_log
        )

        result
      }
    })

    result <- result[!sapply(result, is.null)]
    if (length(result) == 0) {
      return(NULL)
    }

    result
  })

  dist_data <- reactiveVal(list())
  reactivity_paused <- reactiveVal(FALSE)
  byHelpText <- reactiveVal("You can select multiple columns. They will be concatenated if more than one is selected.")
  byHelpTextValue <- reactiveVal("You can select multiple columns. They will be concatenated if more than one is selected.")
  originalByColumns <- reactiveVal(NULL)
  optimalPointResults <- reactiveVal(NULL)
  firstTICalculation <- reactiveVal(FALSE)
  two_sided_ti_initiated <- reactiveVal(FALSE)
  rendering_in_progress <- reactiveVal(list())
  table_rendering <- reactiveVal(FALSE)
  optimalPointResults <- reactiveVal(NULL)
  originalOptimalPointResults <- reactiveVal(NULL)
  single_chain_selected <- reactiveVal(TRUE)
  scatter_checkbox_states <- reactiveVal(list())
  krishnamoorthy_selected_by <- reactiveVal(NULL)


  ### Functions
  inspect_data <- function(data) {
    results <- list()

    tryCatch({
      # Convert to data.frame temporarily to avoid data.table subsetting issues
      is_dt <- is.data.table(data)
      if (is_dt) {
        data_df <- as.data.frame(data)
      } else {
        data_df <- data
      }

      # Check column types
      col_types <- sapply(data_df, class)
      list_cols <- names(which(sapply(col_types, function(x) "list" %in% x)))

      if (length(list_cols) > 0) {
        results$list_warnings <- sprintf("Column '%s' contains list data which may not be suitable for analysis.", list_cols)
      }

      # Check for NA and NaN values
      na_counts <- sapply(data_df, function(col) {
        if (is.list(col)) {
          return(0)  # Skip list columns
        } else {
          return(sum(is.na(col) | is.nan(col)))
        }
      })

      if (sum(na_counts) > 0) {
        results$na_warnings <- sapply(names(na_counts[na_counts > 0]), function(col) {
          sprintf("Column '%s' contains %d NA/NaN values.", col, na_counts[col])
        })
      }

      # Check for columns with both numeric and string data
      mixed_type_cols <- sapply(data_df, function(col) {
        if (is.list(col)) {
          return(FALSE)  # Skip list columns
        } else if (is.character(col)) {
          return(any(grepl("^\\s*-?\\d*\\.?\\d+\\s*$", col)) &&
                   !all(grepl("^\\s*-?\\d*\\.?\\d+\\s*$", col)))
        } else {
          return(FALSE)
        }
      })

      if (any(mixed_type_cols)) {
        results$mixed_type_warnings <- sprintf("Column '%s' contains both numeric and string data.",
                                               names(which(mixed_type_cols)))
      }

      # Check for columns with only integers
      integer_cols <- sapply(data_df, function(col) {
        if (is.list(col) || !is.numeric(col)) {
          return(FALSE)  # Skip list columns and non-numeric columns
        } else {
          return(all(col == floor(col), na.rm = TRUE))
        }
      })

      if (any(integer_cols)) {
        results$integer_warnings <- sprintf("Column '%s' contains only integer values.",
                                            names(which(integer_cols)))
      }

      # Check for columns with only zeros
      zero_cols <- sapply(data_df, function(col) {
        if (is.list(col) || !is.numeric(col)) {
          return(FALSE)  # Skip list columns and non-numeric columns
        } else {
          return(all(col == 0, na.rm = TRUE))
        }
      })

      if (any(zero_cols)) {
        results$zero_warnings <- sprintf("Column '%s' contains only zero values.",
                                         names(which(zero_cols)))
      }

    }, error = function(e) {
      results$warning <- paste("Error during data inspection:", e$message)
    })

    return(results)
  }

  calculate_predicted_distributions <- function() {
    req(InputData(), input$fixed_var, input$random_var, input$by, input$tolerance_level)
    print("Starting predictedDistributions calculation")

    FixedEffects <- input$fixed_var
    RandomParams <- input$random_var

    # Check if "Single Chain Dataset" is selected
    is_single_chain <- identical(input$by, "Single Chain Dataset")
    By <- if(is_single_chain) "SingleChain" else input$by
    ToleranceLevel <- input$tolerance_level

    # If "Single Chain Dataset" is selected, create a SingleChain column if it doesn't exist
    if(By == "SingleChain" && !"SingleChain" %in% names(InputData())) {
      temp_data <- InputData()
      temp_data$SingleChain <- 1
      InputData(temp_data)
    }

    # Create a copy of the data to avoid modifying the original
    dt <- copy(as.data.table(InputData()))

    # Draws must be numeric and complete; NA or infinite draws would otherwise be
    # silently treated as zero. Errors are shown to the user through ErrorMessage.
    check_positive(MultiplicationFactorVariances(), "Multiplication factor")
    check_draw_columns(dt, FixedEffects, "Fixed effect")
    check_draw_columns(dt, RandomParams, "Variance component")
    for (k in 1:3) {
      cov_col_k <- input[[paste0("covariate_", k)]]
      if (!is.null(cov_col_k) && cov_col_k != "None" && cov_col_k %in% names(dt)) {
        check_draw_columns(dt, cov_col_k, "Covariate")
      }
    }
    if (anyNA(dt[[By]])) stop(sprintf("By column '%s' contains missing values", By), call. = FALSE)

    # Print column types for debugging
    print("Column types before processing:")
    for(col in c(FixedEffects, RandomParams)) {
      print(paste(col, ":", class(dt[[col]]), "- Sample values:", paste(head(dt[[col]], 3), collapse=", ")))
    }

    # Initialize mean vector with zeros
    mean_vector <- rep(0, nrow(dt))

    # First, add covariate effects to mean
    for(i in 1:3) {
      covariate_col <- input[[paste0("covariate_", i)]]
      covariate_val <- input[[paste0("covariate_value_", i)]]

      # Only process if a covariate column is selected
      if(!is.null(covariate_col) && covariate_col != "None" && covariate_col %in% names(dt)) {
        # Get the column values
        col_vector <- as.numeric(dt[[covariate_col]])

        # Use covariate value (defaults to 1 if NULL or NA)
        multiplier <- if(is.null(covariate_val) || is.na(covariate_val)) 1 else covariate_val

        # Add covariate effect to mean
        mean_vector <- mean_vector + (col_vector * multiplier)

        print(paste("Added Covariate", i, "effect:", covariate_col, "*", multiplier))
        print(paste("Sample values after adding covariate", i, ":",
                    paste(head(mean_vector, 3), collapse=", ")))
      }
    }

    # Then add fixed effects
    if(length(FixedEffects) > 0) {
      fixed_matrix <- matrix(
        as.numeric(unlist(dt[, ..FixedEffects])),
        nrow = nrow(dt),
        ncol = length(FixedEffects)
      )

      # Add fixed effects to mean
      mean_vector <- mean_vector + rowSums(fixed_matrix, na.rm = TRUE)

      print("Added fixed effects")
      print(paste("Sample values after adding fixed effects:",
                  paste(head(mean_vector, 3), collapse=", ")))
    }

    # Calculate variance
    if(length(RandomParams) > 0) {
      random_matrix <- matrix(
        as.numeric(unlist(dt[, ..RandomParams])),
        nrow = nrow(dt),
        ncol = length(RandomParams)
      )

      # Calculate total variance
      variance_vector <- rowSums(random_matrix, na.rm = TRUE) * MultiplicationFactorVariances()
      check_variance(variance_vector)
    } else {
      variance_vector <- rep(0, nrow(dt))
    }

    # Print matrix types for debugging
    print("Final dimensions and values:")
    print(paste("Number of rows:", length(mean_vector)))
    print(paste("Mean vector sample:", paste(head(mean_vector, 3), collapse=", ")))
    print(paste("Variance vector sample:", paste(head(variance_vector, 3), collapse=", ")))

    # Calculate tolerance parameters
    toleranceaplha <- ToleranceLevel / 100
    toleranceslice <- (1 - toleranceaplha) / 2

    # Calculate quantiles using the mean and variance
    lower_quantiles <- qnorm(toleranceslice, mean = mean_vector, sd = sqrt(pmax(variance_vector, 1e-10)))
    upper_quantiles <- qnorm(toleranceaplha + toleranceslice, mean = mean_vector, sd = sqrt(pmax(variance_vector, 1e-10)))
    medians <- qnorm(0.5, mean = mean_vector, sd = sqrt(pmax(variance_vector, 1e-10)))

    # Create result data.table
    PredictedDistributions_dt <- data.table(
      By = dt[[By]],
      mean = mean_vector,
      median = medians,
      lower_quantile = lower_quantiles,
      upper_quantile = upper_quantiles,
      variance = variance_vector
    )

    # Round to 3 decimals. The inputs are validated to be finite, so every value is
    # a finite number and the columns stay numeric.
    for (col in c("mean", "median", "lower_quantile", "upper_quantile", "variance")) {
      if (any(!is.finite(PredictedDistributions_dt[[col]]))) {
        stop(sprintf("Non-finite values occurred in '%s'; check the scale of the inputs", col),
             call. = FALSE)
      }
      PredictedDistributions_dt[[col]] <- round(PredictedDistributions_dt[[col]], 3)
    }

    # Sort by the "By" column
    setorder(PredictedDistributions_dt, By)

    print("Finished predictedDistributions calculation")
    print(head(PredictedDistributions_dt))

    return(as.data.frame(PredictedDistributions_dt))
  }

  create_distribution_plots <- function(d, id) {
    add_interval_lines <- function(p, d) {
      p %>%
        add_segments(x = d$median, xend = d$median, y = 0, yend = 1, name = "Median", line = list(color = 'black', dash = 'dash'), visible = "legendonly") %>%
        add_segments(x = d$ci_lower, xend = d$ci_lower, y = 0, yend = 1, name = "CI Lower", line = list(color = 'red', dash = 'dot'), visible = "legendonly") %>%
        add_segments(x = d$ci_upper, xend = d$ci_upper, y = 0, yend = 1, name = "CI Upper", line = list(color = 'red', dash = 'dot'), visible = "legendonly") %>%
        add_segments(x = d$pi_lower, xend = d$pi_lower, y = 0, yend = 1, name = "PI Lower", line = list(color = 'green', dash = 'dot'), visible = "legendonly") %>%
        add_segments(x = d$pi_upper, xend = d$pi_upper, y = 0, yend = 1, name = "PI Upper", line = list(color = 'green', dash = 'dot'), visible = "legendonly") %>%
        add_segments(x = d$ti_lower, xend = d$ti_lower, y = 0, yend = 1, name = "TI Lower", line = list(color = 'blue', dash = 'dot'), visible = "legendonly") %>%
        add_segments(x = d$ti_upper, xend = d$ti_upper, y = 0, yend = 1, name = "TI Upper", line = list(color = 'blue', dash = 'dot'), visible = "legendonly")
    }

    normal_plot <- plot_ly(x = d$x, y = d$y_normal, type = 'scatter', mode = 'lines', name = 'Normal') %>%
      layout(title = paste("Normal Distribution for", d$by),
             xaxis = list(title = "Value"),
             yaxis = list(title = "Density")) %>%
      add_interval_lines(d)

    cumulative_plot <- plot_ly(x = d$x, y = d$y_cumulative, type = 'scatter', mode = 'lines', name = 'Cumulative') %>%
      layout(title = paste("Cumulative Distribution for", d$by),
             xaxis = list(title = "Value"),
             yaxis = list(title = "Cumulative Probability")) %>%
      add_interval_lines(d)

    return(list(normal = normal_plot, cumulative = cumulative_plot, id = id))
  }

  calculate_ci <- function(data, ColToAggregate, PercentForCI, ByColumn) {
    # Convert to data.table if not already
    if (!is.data.table(data)) {
      dt <- as.data.table(data)
    } else {
      dt <- data
    }

    if (!ColToAggregate %in% names(dt) || !ByColumn %in% names(dt)) {
      stop(paste("Column", ColToAggregate, "or", ByColumn, "not found in the dataframe"))
    }

    # Calculate alpha and slice once
    alpha <- 1 - PercentForCI/100
    slice <- alpha / 2

    # Use data.table syntax for faster grouping and calculation
    results <- dt[, .(
      N = .N,
      CI_Lower = quantile(get(ColToAggregate), slice),
      CI_Upper = quantile(get(ColToAggregate), 1 - slice),
      Median = median(get(ColToAggregate))
    ), by = ByColumn]

    return(as.data.frame(results))
  }

  calculate_tolerance_interval <- function(data, LowerColToAggregate, UpperColToAggregate, ConfidenceLevel, ByColumn) {
    # Convert to data.table if not already
    if (!is.data.table(data)) {
      dt <- as.data.table(data)
    } else {
      dt <- data
    }

    if (!LowerColToAggregate %in% names(dt) || !UpperColToAggregate %in% names(dt) || !ByColumn %in% names(dt)) {
      stop(paste("One or more required columns not found in the dataframe"))
    }

    # Calculate confidence parameters once
    confidence_alpha <- 1 - ConfidenceLevel/100
    confidence_slice <- confidence_alpha / 2

    # Use data.table syntax for faster grouping and calculation
    results <- dt[, .(
      N = .N,
      Lower_Lower = quantile(get(LowerColToAggregate), confidence_slice),
      Lower_Upper = quantile(get(LowerColToAggregate), 1 - confidence_slice),
      Upper_Lower = quantile(get(UpperColToAggregate), confidence_slice),
      Upper_Upper = quantile(get(UpperColToAggregate), 1 - confidence_slice)
    ), by = ByColumn]

    # Add TI columns
    results[, `:=`(
      TI_Lower = Lower_Lower,
      TI_Upper = Upper_Upper
    )]

    # Select only needed columns
    results <- results[, .(By = get(ByColumn), TI_Lower, TI_Upper)]

    return(as.data.frame(results))
  }

  # Kept in sync with R/calculate_pi.R (see tests/testthat/test-app-parity.R)
  calculate_PI <- function(data, percent_for_pi, Eta = 1e-8) {
    # Check inputs
    if (!is.data.frame(data) || !is.numeric(percent_for_pi) || length(percent_for_pi) != 1 ||
        !is.finite(percent_for_pi) || percent_for_pi <= 0 || percent_for_pi >= 100) {
      stop("Invalid input parameters")
    }
    if (!all(c("By", "mean", "variance") %in% names(data))) {
      stop("'data' must contain the columns By, mean and variance")
    }
    if (!is.numeric(Eta) || length(Eta) != 1 || !is.finite(Eta) || Eta <= 0) {
      stop("'Eta' must be a single positive number")
    }

    # The two probabilities at which the limits are required
    slice <- (1 - percent_for_pi / 100) / 2
    probabilities <- c(slice, 1 - slice)

    # Bisection on  mean_m(Phi((x - mu_m) / sigma_m)) - p
    find_PI_limits <- function(means, variances) {
      sds <- sqrt(pmax(variances, 1e-12))
      vapply(probabilities, function(p) {
        per_draw <- stats::qnorm(p, mean = means, sd = sds)
        lower <- min(per_draw)
        upper <- max(per_draw)
        initial_width <- upper - lower
        # Precision is relative to the spread of the predictive distribution itself
        # (average variance plus variance of the means), not to the size of the bracket
        total_sd <- sqrt(mean(sds^2) + if (length(means) > 1) stats::var(means) else 0)
        tolerance <- Eta * total_sd
        if (initial_width == 0) {
          return(lower)
        }
        # Start from the quantile of the single averaged normal, which is close
        x <- stats::qnorm(p, mean = mean(means), sd = sqrt(mean(sds^2)))
        if (!is.finite(x) || x <= lower || x >= upper) {
          x <- (lower + upper) / 2
        }
        for (iteration in seq_len(200)) {
          # value of (average of the per-draw CDFs) - p, and its bracket update
          value <- mean(stats::pnorm(x, mean = means, sd = sds)) - p
          if (value < 0) {
            lower <- x
          } else {
            upper <- x
          }
          if (upper - lower <= tolerance) {
            break
          }
          # Newton step (the derivative is the average of the per-draw densities);
          # whenever it leaves the bracket, bisect instead, so convergence is guaranteed
          slope <- mean(stats::dnorm(x, mean = means, sd = sds))
          step <- x - value / slope
          if (!is.finite(step) || step <= lower || step >= upper) {
            step <- (lower + upper) / 2
          }
          if (abs(step - x) <= tolerance) {
            x <- step
            break
          }
          x <- step
        }
        x
      }, numeric(1))
    }

    groups <- unique(data$By)
    limits <- lapply(groups, function(group) {
      rows <- which(data$By == group)
      means <- as.numeric(data$mean[rows])
      variances <- as.numeric(data$variance[rows])
      if (anyNA(means) || anyNA(variances)) {
        stop("'mean' and 'variance' must not contain missing values")
      }
      find_PI_limits(means, variances)
    })

    data.frame(
      By = groups,
      PI_Lower = vapply(limits, function(x) x[1], numeric(1)),
      PI_Upper = vapply(limits, function(x) x[2], numeric(1))
    )
  }

  performCalculations <- function(switch_tab = FALSE) {
    tryCatch({
      # Calculate predictedDistributions
      print("Calculating Predicted Distributions")
      incProgress(0.2)
      pred_dist <- calculate_predicted_distributions()
      predictedDistributions(pred_dist)
      print(paste("Predicted Distributions:", nrow(pred_dist), "rows"))

      # Calculate CI
      print("Starting CI calculation")
      incProgress(0.2)
      ci_result <- calculate_ci(pred_dist, "median", input$percent_for_ci, "By")
      print("CI calculation complete")
      print("CI result column names:")
      print(names(ci_result))

      # Calculate TI
      print("Starting TI calculation")
      incProgress(0.2)
      ti_result <- calculate_tolerance_interval(
        pred_dist,
        "lower_quantile",
        "upper_quantile",
        input$confidence_of_tolerance,
        "By"
      )
      print("TI calculation complete")
      print("TI result column names:")
      print(names(ti_result))

      # Calculate PI
      print("Starting PI calculation")
      incProgress(0.2)
      pi_result <- calculate_PI(
        data = pred_dist,
        percent_for_pi = input$percent_for_pi
      )
      print("PI calculation complete")
      print("PI result column names:")
      print(names(pi_result))

      # check column names in PI result for debugging
      print("PI column names:")
      print(names(pi_result))

      # Merge results
      print("Merging results")
      incProgress(0.2)
      combined_result <- merge(merge(ci_result, ti_result, by = "By"), pi_result, by = "By")

      # Failsafe check
      if (any(table(combined_result$By) > 1)) {
        print("Multiple rows detected for some 'By' values. Performing consistency check...")
        check_consistency <- combined_result %>%
          group_by(By) %>%
          summarise(across(everything(), function(x) length(unique(head(x, 5))) == 1)) %>%
          ungroup()

        if (!all(sapply(check_consistency[,-1], all))) {
          warning("Inconsistent values detected for some 'By' groups. Using first row for each group.")
        }

        # Ensure only one row per unique "By" value
        combined_result <- combined_result %>%
          group_by(By) %>%
          summarise(across(everything(), first)) %>%
          ungroup()
      }

      # Add reducible imprecision as percentages (with a notification if it is not
      # defined for a group). With log-normal data it is recalculated afterwards.
      combined_result <- withCallingHandlers(
        add_reducible_imprecision(combined_result, warn = !isTRUE(input$log_normal), percent = TRUE),
        warning = function(w) {
          showNotification(conditionMessage(w), type = "warning")
          invokeRestart("muffleWarning")
        })

      # Store the combined result
      IntervalsOutputInitial(combined_result)
      CICalculated(TRUE)
      print("Calculation complete")

      # Trigger lognorm_adjuster
      lognorm_adjuster()

      # Only switch to the Intervals tab if switch_tab is TRUE
      if (switch_tab) {
        updateTabsetPanel(session, "main_tabs", selected = "One-Sided Intervals")
      }

    }, error = function(e) {
      ErrorMessage(paste("Error in calculating intervals:", e$message))
      print(paste("Error:", e$message))
    })

    if (!is.null(IntervalsOutputActual())) {
      # Calculate plots once and cache them
      interval_plots_cache(render_intervals_plot_one_sided(IntervalsOutputActual(), input))
    }
  }

  lognorm_adjuster <- function() {
    req(IntervalsOutputInitial())
    data <- IntervalsOutputInitial()

    # Handle the Krishnamoorthy TI values if they exist
    if (!is.null(originalOptimalPointResults())) {
      tryCatch({
        # Get the original optimal points
        original_points <- originalOptimalPointResults()

        # Create a modified copy based on log_normal status
        modified_points <- lapply(original_points, function(point) {
          if (!is.null(point$g1) && !is.null(point$g2)) {
            result <- point
            if (input$log_normal) {
              # Apply exp transformation
              result$g1 <- exp(point$g1)
              result$g2 <- exp(point$g2)
            } else {
              # Use original values
              result$g1 <- point$g1
              result$g2 <- point$g2
            }
            return(result)
          } else {
            return(point)
          }
        })

        # Update the reactive value with the modified copy
        optimalPointResults(modified_points)

      }, error = function(e) {
        # Log the error but continue
        print(paste("Error updating optimal points:", e$message))
      })
    }

    # Transform the main data
    if (input$log_normal) {
      # If log_normal is checked, transform the values using exp()
      data <- data %>%
        mutate(across(c(Median, CI_Lower, CI_Upper, TI_Lower, TI_Upper, PI_Lower, PI_Upper), exp))
    }

    # Add reducible imprecision as percentages (notification if not defined)
    data <- withCallingHandlers(
      add_reducible_imprecision(data, percent = TRUE),
      warning = function(w) {
        showNotification(conditionMessage(w), type = "warning")
        invokeRestart("muffleWarning")
      })

    # Update the reactive values
    IntervalsOutputActual(data)
    plot_data(data)

    # Force plot update - this is the key addition
    print("Forcing plot cache update")
    interval_plots_cache(render_intervals_plot_one_sided(IntervalsOutputActual(), input))

    print("lognorm_adjuster completed - plots should update")
  }

  initializeDataTable <- function(data, options = list()) {
    default_options <- list(
      pageLength = 10,
      dom = 't',
      ordering = FALSE,
      searching = FALSE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'))
    )

    combined_options <- modifyList(default_options, options)

    DT::datatable(
      data,
      options = combined_options,
      style = 'bootstrap4',
      rownames = FALSE
    )
  }

  by_updater <- function() {
    by <- input$by
    if (is.null(by) || length(by) == 0) {
      return(list(selected = "Single Chain Dataset", choices = c("Single Chain Dataset", names(InputData())), original_columns = NULL))
    } else if (length(by) > 1) {
      # Merge the selected columns
      data <- InputData()
      merged_column <- apply(data[, by], 1, paste, collapse = "_")

      # Create new column name
      new_column_name <- "Merged_By_Columns"

      # Add the new column to InputData
      data[[new_column_name]] <- merged_column
      InputData(data)

      # Print first 5 values of the new column
      print(paste("First 5 values of", new_column_name, ":"))
      print(head(merged_column, 5))

      # Return the new column name, updated choices, and original columns
      return(list(selected = new_column_name, choices = c("Single Chain Dataset", names(data)), original_columns = by))
    } else {
      return(list(selected = by[1], choices = c("Single Chain Dataset", names(InputData())), original_columns = NULL))
    }
  }

  parse_x_value_input <- function(input_text) {
    # Check if input is null or NA
    if (is.null(input_text) || is.na(input_text)) {
      return(NULL)
    }

    # Check if input matches range pattern [value1,value2]
    # This regex handles spaces around the comma and values
    range_pattern <- "^\\s*\\[\\s*(-?\\d*\\.?\\d+)\\s*,\\s*(-?\\d*\\.?\\d+)\\s*\\]\\s*$"

    if (grepl(range_pattern, input_text)) {
      # Extract the two values from the range notation
      matches <- regmatches(input_text, regexec(range_pattern, input_text))[[1]]
      xmin <- as.numeric(matches[2])
      xmax <- as.numeric(matches[3])

      # Ensure xmin is less than xmax
      if (xmin > xmax) {
        temp <- xmin
        xmin <- xmax
        xmax <- temp
      }

      return(list(type = "range", xmin = xmin, xmax = xmax))
    } else {
      # Try to parse as a single numeric value
      tryCatch({
        value <- as.numeric(input_text)
        if (!is.na(value)) {
          return(list(type = "single", value = value))
        } else {
          return(NULL)
        }
      }, error = function(e) {
        return(NULL)
      })
    }
  }

  render_intervals_table <- function(ci_data, input, is_one_sided = FALSE) {
    print("Rendering intervals table")
    ci_percent <- input$percent_for_ci
    pi_percent <- input$percent_for_pi
    ti_percent <- input$tolerance_level
    ti_confidence <- input$confidence_of_tolerance

    is_single_chain <- single_chain_selected() && length(unique(ci_data$By)) == 1

    # Get the actual column name used for "by"
    by_column_name <- if(input$by == "Single Chain Dataset") {
      "By Identifier"
    } else if(input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
      paste(originalByColumns(), collapse = ", ")
    } else {
      input$by
    }

    if (nrow(ci_data) == 0) {
      return(datatable(data.frame(Message = "No data available for the selected criteria."),
                       options = list(dom = 't'),
                       style = 'bootstrap4'))
    }

    # Check if columns exist and standardize names if needed
    if("PI_lower" %in% names(ci_data) && !"PI_Lower" %in% names(ci_data)) {
      ci_data$PI_Lower <- ci_data$PI_lower
      ci_data$PI_Upper <- ci_data$PI_upper
    }

    # Create a display version of the data for single chain case
    display_data <- ci_data
    if (is_single_chain) {
      display_data$By <- "" # Empty the By column for display purposes
    }

    # Define columns and headers based on whether it's one-sided or two-sided
    if(is_one_sided) {
      # ONE-SIDED TABLE (keep all columns including PI)
      if(input$show_reducible){
        display_columns <- c("By", "Median", "CI_Lower", "CI_Upper", "PI_Lower", "PI_Upper", "TI_Lower", "TI_Upper", "ReducibleLower", "ReducibleUpper")

        # Check if all columns exist
        missing_cols <- setdiff(display_columns, names(ci_data))
        if(length(missing_cols) > 0) {
          print(paste("Missing columns:", paste(missing_cols, collapse=", ")))
          return(datatable(data.frame(Message = paste("Missing columns in data:", paste(missing_cols, collapse=", "))),
                           options = list(dom = 't'),
                           style = 'bootstrap4'))
        }

        column_names <- if (is_single_chain) {
          c("", "Median", "CI Lower", "CI Upper", "PI Lower", "PI Upper", "TI Lower", "TI Upper", "Reducible Lower", "Reducible Upper")
        } else {
          c(by_column_name, "Median", "CI Lower", "CI Upper", "PI Lower", "PI Upper", "TI Lower", "TI Upper", "Reducible Lower", "Reducible Upper")
        }

        header_names <- c(
          "By Identifier",
          "Point Estimate",
          sprintf("One-Sided %s%% CI", ci_percent),
          sprintf("One-Sided %s%% PI", pi_percent),
          sprintf("One-Sided %s%% TI (%s%% Credibility)", ti_percent, ti_confidence),
          "Reducible Imprecision (%)"
        )

        header_colspans <- c(1, 1, 2, 2, 2, 2)

        container <- tags$table(
          class = 'display',
          tags$thead(
            tags$tr(
              tags$th(rowspan = 2, if (is_single_chain) "" else by_column_name),
              tags$th(rowspan = 2, "Point Estimate"),
              do.call(
                tagList,
                mapply(function(name, colspan) tags$th(colspan = colspan, name),
                       header_names[-c(1, 2)], header_colspans[-c(1, 2)], SIMPLIFY = FALSE)
              )
            ),
            tags$tr(
              tags$th("Lower"), tags$th("Upper"), # CI
              tags$th("Lower"), tags$th("Upper"), # PI
              tags$th("Lower"), tags$th("Upper"), # TI
              tags$th("Lower"), tags$th("Upper")  # Reducible
            )
          )
        )
      } else {
        display_columns <- c("By", "Median", "CI_Lower", "CI_Upper", "PI_Lower", "PI_Upper", "TI_Lower", "TI_Upper")

        # Check if all columns exist
        missing_cols <- setdiff(display_columns, names(ci_data))
        if(length(missing_cols) > 0) {
          print(paste("Missing columns:", paste(missing_cols, collapse=", ")))
          return(datatable(data.frame(Message = paste("Missing columns in data:", paste(missing_cols, collapse=", "))),
                           options = list(dom = 't'),
                           style = 'bootstrap4'))
        }

        column_names <- if (is_single_chain) {
          c("", "Median", "CI Lower", "CI Upper", "PI Lower", "PI Upper", "TI Lower", "TI Upper")
        } else {
          c(by_column_name, "Median", "CI Lower", "CI Upper", "PI Lower", "PI Upper", "TI Lower", "TI Upper")
        }

        header_names <- c(
          "By Identifier",
          "Point Estimate",
          sprintf("One-Sided %s%% CI", ci_percent),
          sprintf("One-Sided %s%% PI", pi_percent),
          sprintf("One-Sided %s%% TI (%s%% Credibility)", ti_percent, ti_confidence)
        )

        header_colspans <- c(1, 1, 2, 2, 2)

        container <- tags$table(
          class = 'display',
          tags$thead(
            tags$tr(
              tags$th(rowspan = 2, if (is_single_chain) "" else by_column_name),
              tags$th(rowspan = 2, "Point Estimate"),
              do.call(
                tagList,
                mapply(function(name, colspan) tags$th(colspan = colspan, name),
                       header_names[-c(1, 2)], header_colspans[-c(1, 2)], SIMPLIFY = FALSE)
              )
            ),
            tags$tr(
              tags$th("Lower"), tags$th("Upper"), # CI
              tags$th("Lower"), tags$th("Upper"), # PI
              tags$th("Lower"), tags$th("Upper")  # TI
            )
          )
        )
      }
    } else {
      # TWO-SIDED TABLE (exclude PI columns)
      display_columns <- c("By", "Median", "CI_Lower", "CI_Upper", "TI_Lower", "TI_Upper")

      # Check if all columns exist
      missing_cols <- setdiff(display_columns, names(ci_data))
      if(length(missing_cols) > 0) {
        print(paste("Missing columns:", paste(missing_cols, collapse=", ")))
        return(datatable(data.frame(Message = paste("Missing columns in data:", paste(missing_cols, collapse=", "))),
                         options = list(dom = 't'),
                         style = 'bootstrap4'))
      }

      column_names <- if (is_single_chain) {
        c("", "Median", "CI Lower", "CI Upper", "TI Lower", "TI Upper")
      } else {
        c(by_column_name, "Median", "CI Lower", "CI Upper", "TI Lower", "TI Upper")
      }

      header_names <- c(
        "By Identifier",
        "Point Estimate",
        sprintf("Two-Sided %s%% CI", ci_percent),
        sprintf("Two-Sided %s%% TI (%s%% Credibility)", ti_percent, ti_confidence)
      )

      header_colspans <- c(1, 1, 2, 2)

      container <- tags$table(
        class = 'display',
        tags$thead(
          tags$tr(
            tags$th(rowspan = 2, if (is_single_chain) "" else by_column_name),
            tags$th(rowspan = 2, "Point Estimate"),
            do.call(
              tagList,
              mapply(function(name, colspan) tags$th(colspan = colspan, name),
                     header_names[-c(1, 2)], header_colspans[-c(1, 2)], SIMPLIFY = FALSE)
            )
          ),
          tags$tr(
            tags$th("Lower"), tags$th("Upper"), # CI
            tags$th("Lower"), tags$th("Upper")  # TI
          )
        )
      )
    }

    # Add CSS for padding directly to the output options
    datatable(display_data[, display_columns],
              options = list(
                pageLength = 50,
                dom = 't',
                scrollX = TRUE,
                scrollY = "400px",
                fixedHeader = TRUE,
                scrollCollapse = TRUE,
                paging = FALSE,
                ordering = TRUE,
                order = list(),
                columnDefs = list(
                  list(className = 'dt-center', targets = '_all'),
                  list(className = 'dt-left', targets = 0)
                )
              ),
              style = 'bootstrap4',
              container = container,
              rownames = FALSE
    ) %>%
      formatRound(columns = setdiff(display_columns, "By"), digits = 4)
  }

  render_intervals_plot <- function(data, input) {
    print("Rendering intervals plot")

    tryCatch({
      if(nrow(data) == 0) {
        stop("No data available for plotting")
      }

      # Check for column names and standardize if needed
      if("PI_lower" %in% names(data) && !"PI_Lower" %in% names(data)) {
        data$PI_Lower <- data$PI_lower
        data$PI_Upper <- data$PI_upper
      }

      # Check that all required columns exist
      required_cols <- c("By", "CI_Lower", "CI_Upper", "TI_Lower", "TI_Upper", "PI_Lower", "PI_Upper", "Median")
      missing_cols <- setdiff(required_cols, names(data))
      if(length(missing_cols) > 0) {
        stop(paste("Missing required columns in the data:", paste(missing_cols, collapse=", ")))
      }

      # Convert By to factor and get number of unique values
      data$By <- factor(data$By)
      n_by <- length(unique(data$By))

      # Check if this is a single chain situation
      is_single_chain <- single_chain_selected() && n_by == 1 && (as.character(data$By[1]) == "1" || as.character(data$By[1]) == "SingleChain")

      # Determine x-axis title based on the By values
      x_axis_title <- if (is_single_chain) {
        "By Identifier"
      } else if (n_by == 1) {
        # If there's only one By value and it's not SingleChain, use that value
        as.character(data$By[1])
      } else if (input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
        # If using merged columns, show the original column names
        paste(originalByColumns(), collapse = ", ")
      } else if (length(input$by) == 1 && input$by != "Single Chain Dataset") {
        # If a single column is selected, use that column name
        input$by
      } else {
        # Default case
        "By Identifier"
      }

      # Improved scaling logic for cap width and offsets
      # Base values that work well for medium number of points
      base_cap_width <- 0.05
      base_offset <- 0.18  # Increased from 0.15 to create more space between intervals

      # Calculate adaptive scaling factors with better constraints
      # For cap width: starts larger for few points, gets smaller but never too small
      cap_scaling <- 5 / (n_by + 4)  # Smoother decline curve
      cap_width <- base_cap_width * cap_scaling
      cap_width <- max(0.02, min(0.1, cap_width))  # Constrain between 0.02 and 0.1

      # For offset: similar adaptive scaling but with different constraints
      offset_scaling <- 4 / (n_by + 3)  # Smoother decline for offsets
      offset_scale <- base_offset * offset_scaling
      offset_scale <- max(0.08, min(0.2, offset_scale))  # Increased minimum from 0.05 to 0.08

      # Initialize the plot with better spacing for x-axis
      p <- plot_ly()

      # Add reducible uncertainty ribbons first (if enabled) so they appear behind the lines
      if (input$show_reducible) {
        # Create a data frame for the ribbons
        ribbon_data <- data.frame()

        for (i in 1:nrow(data)) {
          by_val <- as.numeric(data$By[i])
          by_label <- as.character(data$By[i])

          # Width of the ribbon (centered on the PI position)
          ribbon_width <- 0.15

          # Add rows for upper reducible uncertainty
          upper_ribbon <- data.frame(
            x = c(by_val - ribbon_width, by_val + ribbon_width),
            y = data$PI_Upper[i],
            ymin = data$PI_Upper[i],
            ymax = data$TI_Upper[i],
            by_val = by_val,
            by_label = by_label,
            type = "Upper",
            reducible_pct = data$ReducibleUpper[i]
          )

          # Add rows for lower reducible uncertainty
          lower_ribbon <- data.frame(
            x = c(by_val - ribbon_width, by_val + ribbon_width),
            y = data$PI_Lower[i],
            ymin = data$TI_Lower[i],
            ymax = data$PI_Lower[i],
            by_val = by_val,
            by_label = by_label,
            type = "Lower",
            reducible_pct = data$ReducibleLower[i]
          )

          ribbon_data <- rbind(ribbon_data, upper_ribbon, lower_ribbon)
        }

        # Add upper reducible uncertainty ribbons
        upper_ribbons <- ribbon_data[ribbon_data$type == "Upper",]
        if (nrow(upper_ribbons) > 0) {
          for (i in 1:nrow(upper_ribbons)) {
            p <- add_ribbons(p,
                             x = c(upper_ribbons$x[2*i-1], upper_ribbons$x[2*i]),
                             ymin = rep(upper_ribbons$ymin[2*i-1], 2),
                             ymax = rep(upper_ribbons$ymax[2*i-1], 2),
                             fillcolor = 'rgba(255, 165, 0, 0.6)',  # Changed to opaque orange
                             line = list(color = 'rgba(255, 165, 0, 1)', width = 1),  # Fully opaque line
                             name = "Reducible Imprecision",  # Consistent name
                             showlegend = (i == 1),  # Only show first one in legend
                             hoverinfo = "text",
                             text = sprintf("Reducible Imprecision<br>By: %s<br>Value: %.2f%%<br>Range: %.4f to %.4f",
                                            upper_ribbons$by_label[2*i-1],
                                            upper_ribbons$reducible_pct[2*i-1],
                                            upper_ribbons$ymin[2*i-1],
                                            upper_ribbons$ymax[2*i-1]))
          }
        }

        # Add lower reducible uncertainty ribbons
        lower_ribbons <- ribbon_data[ribbon_data$type == "Lower",]
        if (nrow(lower_ribbons) > 0) {
          for (i in 1:nrow(lower_ribbons)) {
            p <- add_ribbons(p,
                             x = c(lower_ribbons$x[2*i-1], lower_ribbons$x[2*i]),
                             ymin = rep(lower_ribbons$ymin[2*i-1], 2),
                             ymax = rep(lower_ribbons$ymax[2*i-1], 2),
                             fillcolor = 'rgba(255, 165, 0, 0.6)',  # Changed to opaque orange
                             line = list(color = 'rgba(255, 165, 0, 1)', width = 1),  # Fully opaque line
                             name = "Reducible Imprecision",  # Consistent name
                             showlegend = FALSE,  # Never show in legend
                             hoverinfo = "text",
                             text = sprintf("Reducible Imprecision<br>By: %s<br>Value: %.2f%%<br>Range: %.4f to %.4f",
                                            lower_ribbons$by_label[2*i-1],
                                            lower_ribbons$reducible_pct[2*i-1],
                                            lower_ribbons$ymin[2*i-1],
                                            lower_ribbons$ymax[2*i-1]))
          }
        }
      }

      # Add traces for each interval type and By value
      for (i in 1:nrow(data)) {
        by_val <- as.numeric(data$By[i])
        by_label <- as.character(data$By[i])

        # Add CI intervals
        p <- add_segments(p,
                          x = by_val - offset_scale, y = data$CI_Lower[i],
                          xend = by_val - offset_scale, yend = data$CI_Upper[i],
                          line = list(color = '#de2f2f', width = 1),
                          showlegend = i == 1,
                          name = "CI",
                          hoverinfo = "text",
                          text = sprintf("CI<br>By: %s<br>Lower: %.4f<br>Upper: %.4f<br>Median: %.4f",
                                         by_label, data$CI_Lower[i], data$CI_Upper[i], data$Median[i]))

        # Add CI caps
        p <- add_segments(p,
                          x = by_val - offset_scale - cap_width, y = data$CI_Lower[i],
                          xend = by_val - offset_scale + cap_width, yend = data$CI_Lower[i],
                          line = list(color = '#de2f2f', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("CI Lower: %.4f", data$CI_Lower[i]))

        p <- add_segments(p,
                          x = by_val - offset_scale - cap_width, y = data$CI_Upper[i],
                          xend = by_val - offset_scale + cap_width, yend = data$CI_Upper[i],
                          line = list(color = '#de2f2f', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("CI Upper: %.4f", data$CI_Upper[i]))

        # Add CI median point
        p <- add_markers(p,
                         x = by_val - offset_scale, y = data$Median[i],
                         marker = list(color = '#de2f2f', size = 3.5),
                         showlegend = FALSE,
                         hoverinfo = "text",
                         text = sprintf("CI Median: %.4f", data$Median[i]))

        # Add PI intervals
        p <- add_segments(p,
                          x = by_val, y = data$PI_Lower[i],
                          xend = by_val, yend = data$PI_Upper[i],
                          line = list(color = '#006400', width = 1),
                          showlegend = i == 1,
                          name = "PI",
                          hoverinfo = "text",
                          text = sprintf("PI<br>By: %s<br>Lower: %.4f<br>Upper: %.4f<br>Median: %.4f",
                                         by_label, data$PI_Lower[i], data$PI_Upper[i], data$Median[i]))

        # Add PI caps
        p <- add_segments(p,
                          x = by_val - cap_width, y = data$PI_Lower[i],
                          xend = by_val + cap_width, yend = data$PI_Lower[i],
                          line = list(color = '#006400', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("PI Lower: %.4f", data$PI_Lower[i]))

        p <- add_segments(p,
                          x = by_val - cap_width, y = data$PI_Upper[i],
                          xend = by_val + cap_width, yend = data$PI_Upper[i],
                          line = list(color = '#006400', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("PI Upper: %.4f", data$PI_Upper[i]))

        # Add PI median point
        p <- add_markers(p,
                         x = by_val, y = data$Median[i],
                         marker = list(color = '#006400', size = 3.5),
                         showlegend = FALSE,
                         hoverinfo = "text",
                         text = sprintf("PI Median: %.4f", data$Median[i]))

        # Add TI intervals
        p <- add_segments(p,
                          x = by_val + offset_scale, y = data$TI_Lower[i],
                          xend = by_val + offset_scale, yend = data$TI_Upper[i],
                          line = list(color = '#4232e7', width = 1),
                          showlegend = i == 1,
                          name = "TI",
                          hoverinfo = "text",
                          text = sprintf("TI<br>By: %s<br>Lower: %.4f<br>Upper: %.4f<br>Median: %.4f",
                                         by_label, data$TI_Lower[i], data$TI_Upper[i], data$Median[i]))

        # Add TI caps
        p <- add_segments(p,
                          x = by_val + offset_scale - cap_width, y = data$TI_Lower[i],
                          xend = by_val + offset_scale + cap_width, yend = data$TI_Lower[i],
                          line = list(color = '#4232e7', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("TI Lower: %.4f", data$TI_Lower[i]))

        p <- add_segments(p,
                          x = by_val + offset_scale - cap_width, y = data$TI_Upper[i],
                          xend = by_val + offset_scale + cap_width, yend = data$TI_Upper[i],
                          line = list(color = '#4232e7', width = 1.5),
                          showlegend = FALSE,
                          hoverinfo = "text",
                          text = sprintf("TI Upper: %.4f", data$TI_Upper[i]))

        # Add TI median point
        p <- add_markers(p,
                         x = by_val + offset_scale, y = data$Median[i],
                         marker = list(color = '#4232e7', size = 3.5),
                         showlegend = FALSE,
                         hoverinfo = "text",
                         text = sprintf("TI Median: %.4f", data$Median[i]))
      }

      # Improved x-axis layout with better spacing
      # Calculate range padding based on number of points - now applied to both sides
      range_padding <- max(0.8, min(2.5, 6/n_by))

      # For single chain, replace the "1" with an empty string
      tick_text <- levels(data$By)
      if (is_single_chain) {
        tick_text <- ""  # Replace with empty string for single chain
      }

      p <- layout(p,
                  title = "Bayesian Statistical Intervals",  # Updated title
                  xaxis = list(
                    title = x_axis_title,  # Dynamic x-axis title
                    tickmode = "array",
                    tickvals = 1:length(levels(data$By)),
                    ticktext = tick_text,
                    tickangle = 45,
                    range = c(1 - range_padding, length(levels(data$By)) + range_padding)  # Dynamic range with padding on both sides
                  ),
                  yaxis = list(title = "Value"),
                  legend = list(x = 1.05, y = 0.5),
                  margin = list(b = 100)  # Add bottom margin for rotated labels
      )

      return(p)

    }, error = function(e) {
      print(paste("Error in plotting:", e$message))
      plotly_empty() %>%
        layout(title = paste("Error:", e$message))
    })
  }

  render_intervals_plot_one_sided <- function(data, input) {
    print("Rendering one-sided intervals plot")

    tryCatch({
      if(nrow(data) == 0) {
        stop("No data available for plotting")
      }

      # Check for column names and standardize if needed
      if("PI_lower" %in% names(data) && !"PI_Lower" %in% names(data)) {
        data$PI_Lower <- data$PI_lower
        data$PI_Upper <- data$PI_upper
      }

      # Check that all required columns exist
      required_cols <- c("By", "CI_Lower", "CI_Upper", "TI_Lower", "TI_Upper", "PI_Lower", "PI_Upper", "Median")
      missing_cols <- setdiff(required_cols, names(data))
      if(length(missing_cols) > 0) {
        stop(paste("Missing required columns in the data:", paste(missing_cols, collapse=", ")))
      }

      # Convert By to factor and get number of unique values
      data$By <- factor(data$By)
      n_by <- length(unique(data$By))

      # Check if this is a single chain situation
      is_single_chain <- single_chain_selected() && n_by == 1 && (as.character(data$By[1]) == "1" || as.character(data$By[1]) == "SingleChain")

      # Determine x-axis title based on the By values
      x_axis_title <- if (is_single_chain) {
        "By Identifier"
      } else if (n_by == 1) {
        as.character(data$By[1])
      } else if (input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
        paste(originalByColumns(), collapse = ", ")
      } else if (length(input$by) == 1 && input$by != "Single Chain Dataset") {
        input$by
      } else {
        "By Identifier"
      }

      # Calculate plot range for extending to infinity
      all_values <- c(data$CI_Lower, data$CI_Upper, data$PI_Lower, data$PI_Upper, data$TI_Lower, data$TI_Upper)
      y_range <- range(all_values, na.rm = TRUE)
      y_extend <- diff(y_range) * 0.3  # Extend by 30% of the range
      y_min_extended <- y_range[1] - y_extend
      y_max_extended <- y_range[2] + y_extend

      # Improved scaling logic for cap width and offsets
      base_cap_width <- 0.05
      base_offset <- 0.18

      cap_scaling <- 5 / (n_by + 4)
      cap_width <- base_cap_width * cap_scaling
      cap_width <- max(0.02, min(0.1, cap_width))

      offset_scaling <- 4 / (n_by + 3)
      offset_scale <- base_offset * offset_scaling
      offset_scale <- max(0.08, min(0.2, offset_scale))

      # Create separate plots for upper and lower bounds
      upper_plot <- plot_ly() %>% layout(title = "Upper Bounds")
      lower_plot <- plot_ly() %>% layout(title = "Lower Bounds")

      # Add traces for upper bounds plot
      for (i in 1:nrow(data)) {
        by_val <- as.numeric(data$By[i])
        by_label <- as.character(data$By[i])

        # Add reducible uncertainty ribbons if enabled
        if (input$show_reducible) {
          ribbon_width <- 0.15
          # Upper plot ribbons
          upper_plot <- upper_plot %>%
            add_ribbons(
              x = c(by_val - ribbon_width, by_val + ribbon_width),
              ymin = rep(data$PI_Upper[i], 2),
              ymax = rep(data$TI_Upper[i], 2),
              fillcolor = 'rgba(255, 165, 0, 0.6)',
              line = list(color = 'rgba(255, 165, 0, 1)', width = 1),
              name = "Reducible Imprecision",
              showlegend = (i == 1),
              hoverinfo = "text",
              text = sprintf("Reducible Imprecision<br>By: %s<br>Range: %.4f to %.4f",
                             by_label, data$PI_Upper[i], data$TI_Upper[i])
            )

          # Lower plot ribbons
          lower_plot <- lower_plot %>%
            add_ribbons(
              x = c(by_val - ribbon_width, by_val + ribbon_width),
              ymin = rep(data$TI_Lower[i], 2),
              ymax = rep(data$PI_Lower[i], 2),
              fillcolor = 'rgba(255, 165, 0, 0.6)',
              line = list(color = 'rgba(255, 165, 0, 1)', width = 1),
              name = "Reducible Imprecision",
              showlegend = (i == 1),
              hoverinfo = "text",
              text = sprintf("Reducible Imprecision<br>By: %s<br>Range: %.4f to %.4f",
                             by_label, data$TI_Lower[i], data$PI_Lower[i])
            )
        }

        # Add CI Upper traces
        upper_plot <- upper_plot %>%
          add_segments(
            x = by_val - offset_scale, y = data$Median[i],
            xend = by_val - offset_scale, yend = data$CI_Upper[i],
            line = list(color = '#de2f2f', width = 2),
            showlegend = i == 1,
            name = "CI Upper"
          ) %>%
          add_segments(
            x = by_val - offset_scale, y = data$CI_Upper[i],
            xend = by_val - offset_scale, yend = y_min_extended,
            line = list(color = '#de2f2f', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val - offset_scale - cap_width, y = data$CI_Upper[i],
            xend = by_val - offset_scale + cap_width, yend = data$CI_Upper[i],
            line = list(color = '#de2f2f', width = 2),
            showlegend = FALSE
          )

        # Add PI Upper traces
        upper_plot <- upper_plot %>%
          add_segments(
            x = by_val, y = data$Median[i],
            xend = by_val, yend = data$PI_Upper[i],
            line = list(color = '#006400', width = 2),
            showlegend = i == 1,
            name = "PI Upper"
          ) %>%
          add_segments(
            x = by_val, y = data$PI_Upper[i],
            xend = by_val, yend = y_min_extended,
            line = list(color = '#006400', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val - cap_width, y = data$PI_Upper[i],
            xend = by_val + cap_width, yend = data$PI_Upper[i],
            line = list(color = '#006400', width = 2),
            showlegend = FALSE
          )

        # Add TI Upper traces
        upper_plot <- upper_plot %>%
          add_segments(
            x = by_val + offset_scale, y = data$Median[i],
            xend = by_val + offset_scale, yend = data$TI_Upper[i],
            line = list(color = '#4232e7', width = 2),
            showlegend = i == 1,
            name = "TI Upper"
          ) %>%
          add_segments(
            x = by_val + offset_scale, y = data$TI_Upper[i],
            xend = by_val + offset_scale, yend = y_min_extended,
            line = list(color = '#4232e7', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val + offset_scale - cap_width, y = data$TI_Upper[i],
            xend = by_val + offset_scale + cap_width, yend = data$TI_Upper[i],
            line = list(color = '#4232e7', width = 2),
            showlegend = FALSE
          )

        # Add median points for upper plot
        upper_plot <- upper_plot %>%
          add_markers(
            x = c(by_val - offset_scale, by_val, by_val + offset_scale),
            y = rep(data$Median[i], 3),
            marker = list(color = c('#de2f2f', '#006400', '#4232e7'), size = 4),
            showlegend = FALSE
          )

        # Add traces for lower bounds plot with correct data points
        lower_plot <- lower_plot %>%
          # CI Lower traces
          add_segments(
            x = by_val - offset_scale, y = data$Median[i],
            xend = by_val - offset_scale, yend = data$CI_Lower[i],
            line = list(color = '#de2f2f', width = 2),
            showlegend = i == 1,
            name = "CI Lower"
          ) %>%
          add_segments(
            x = by_val - offset_scale, y = data$CI_Lower[i],
            xend = by_val - offset_scale, yend = y_max_extended,
            line = list(color = '#de2f2f', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val - offset_scale - cap_width, y = data$CI_Lower[i],
            xend = by_val - offset_scale + cap_width, yend = data$CI_Lower[i],
            line = list(color = '#de2f2f', width = 2),
            showlegend = FALSE
          ) %>%

          # PI Lower traces
          add_segments(
            x = by_val, y = data$Median[i],
            xend = by_val, yend = data$PI_Lower[i],
            line = list(color = '#006400', width = 2),
            showlegend = i == 1,
            name = "PI Lower"
          ) %>%
          add_segments(
            x = by_val, y = data$PI_Lower[i],
            xend = by_val, yend = y_max_extended,
            line = list(color = '#006400', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val - cap_width, y = data$PI_Lower[i],
            xend = by_val + cap_width, yend = data$PI_Lower[i],
            line = list(color = '#006400', width = 2),
            showlegend = FALSE
          ) %>%

          # TI Lower traces
          add_segments(
            x = by_val + offset_scale, y = data$Median[i],
            xend = by_val + offset_scale, yend = data$TI_Lower[i],
            line = list(color = '#4232e7', width = 2),
            showlegend = i == 1,
            name = "TI Lower"
          ) %>%
          add_segments(
            x = by_val + offset_scale, y = data$TI_Lower[i],
            xend = by_val + offset_scale, yend = y_max_extended,
            line = list(color = '#4232e7', width = 2),
            showlegend = FALSE
          ) %>%
          add_segments(
            x = by_val + offset_scale - cap_width, y = data$TI_Lower[i],
            xend = by_val + offset_scale + cap_width, yend = data$TI_Lower[i],
            line = list(color = '#4232e7', width = 2),
            showlegend = FALSE
          ) %>%

          # Add median points for lower plot
          add_markers(
            x = c(by_val - offset_scale, by_val, by_val + offset_scale),
            y = rep(data$Median[i], 3),
            marker = list(color = c('#de2f2f', '#006400', '#4232e7'), size = 4),
            showlegend = FALSE
          )
      }

      # Calculate range padding
      range_padding <- max(0.8, min(2.5, 6/n_by))

      # For single chain, replace the "1" with an empty string
      tick_text <- levels(data$By)
      if (is_single_chain) {
        tick_text <- ""
      }

      # Set layout for both plots
      upper_plot <- upper_plot %>% layout(
        xaxis = list(
          title = x_axis_title,
          tickmode = "array",
          tickvals = 1:length(levels(data$By)),
          ticktext = tick_text,
          tickangle = 45,
          range = c(1 - range_padding, length(levels(data$By)) + range_padding)
        ),
        yaxis = list(title = "Upper Bound Values"),
        legend = list(x = 1.05, y = 0.5),
        margin = list(b = 100)
      )

      lower_plot <- lower_plot %>% layout(
        xaxis = list(
          title = x_axis_title,
          tickmode = "array",
          tickvals = 1:length(levels(data$By)),
          ticktext = tick_text,
          tickangle = 45,
          range = c(1 - range_padding, length(levels(data$By)) + range_padding)
        ),
        yaxis = list(title = "Lower Bound Values"),
        legend = list(x = 1.05, y = 0.5),
        margin = list(b = 100)
      )

      return(list(upper = upper_plot, lower = lower_plot))

    }, error = function(e) {
      print(paste("Error in plotting:", e$message))
      return(list(
        upper = plotly_empty() %>% layout(title = paste("Error:", e$message)),
        lower = plotly_empty() %>% layout(title = paste("Error:", e$message))
      ))
    })
  }

  calculate_gP_mean <- function(scatter_data, P) {
    # Ensure scatter_data is a data frame
    scatter_data <- as.data.frame(scatter_data)

    # Check if required columns exist
    required_cols <- c("By", "lower_quantile", "upper_quantile")
    if (!all(required_cols %in% names(scatter_data))) {
      stop("Required columns (By, lower_quantile, upper_quantile) not found in the data")
    }

    # Calculate gP_mean for each By value
    gP_mean <- scatter_data %>%
      group_by(By) %>%
      summarise(
        M = n(),  # number of samples
        gP_mean = (1/M) * sum((lower_quantile + upper_quantile)/2),
        .groups = 'drop'
      )

    # Convert to a named vector for easier access
    gP_mean_vector <- setNames(gP_mean$gP_mean, gP_mean$By)

    return(gP_mean_vector)
  }

  calculate_gamma_difference <- function(scatter_data, gP_mean, point_x_vector, gamma_actual) {
    # Use the current confidence level from input
    gamma_actual <- isolate(input$confidence_of_tolerance)

    # Vectorized calculations
    point_y_vector <- 2 * gP_mean - point_x_vector
    total_points <- nrow(scatter_data)

    # Create a matrix of conditions
    conditions <- outer(scatter_data$upper_quantile, point_y_vector, "<=") &
      outer(scatter_data$lower_quantile, point_x_vector, ">=")

    # Sum the conditions for each point_x
    points_satisfying <- colSums(conditions)

    gamma_for_points <- (points_satisfying / total_points) * 100
    gamma_differences <- gamma_actual - gamma_for_points

    return(list(
      gamma_differences = gamma_differences,
      points_satisfying = points_satisfying,
      total_points = total_points
    ))
  }

  find_optimal_point <- function(scatter_data, gP_mean_vector, gamma_actual, max_iterations = 300, tolerance = 0.01) {
    # Initialize results list
    results <- list()

    # Process each By value separately
    by_values <- unique(scatter_data$By)

    for(by_val in by_values) {
      cat("\n\n=== Processing By Value:", by_val, "===\n")

      # Get data for this By value
      by_data <- scatter_data[scatter_data$By == by_val, ]

      # Get gP_mean for this By value
      by_gP_mean <- gP_mean_vector[as.character(by_val)]
      cat("gP_mean for this group:", by_gP_mean, "\n")

      # Get min and max x values for this By value
      x_min <- min(by_data$lower_quantile)
      x_max <- max(by_data$lower_quantile)
      cat("Search range: [", x_min, ",", x_max, "]\n")

      # Initialize bisection variables
      a <- x_min
      b <- x_max
      iteration <- 0

      # Bisection method
      while(iteration < max_iterations) {
        # Calculate midpoint
        c <- (a + b) / 2

        # Calculate gamma differences at a and c
        diff_a <- calculate_gamma_difference(by_data, by_gP_mean, a, gamma_actual)
        diff_c <- calculate_gamma_difference(by_data, by_gP_mean, c, gamma_actual)

        if(iteration %% 10 == 0) { # Print every 10 iterations
          cat("\nIteration", iteration, ":\n")
          cat("Current interval: [", a, ",", b, "]\n")
          cat("Midpoint:", c, "\n")
          cat("Gamma at midpoint:", sprintf("%.2f", 100 - diff_c$gamma_difference), "% (Target:", gamma_actual, "%)\n")
          cat("Difference from target:", sprintf("%.2f", diff_c$gamma_difference), "\n")
        }

        # Check if we're close enough
        if(abs(diff_c$gamma_difference) < tolerance) {
          cat("\nConverged! Difference from target (", sprintf("%.2f", diff_c$gamma_difference),
              ") is within tolerance (", tolerance, ")\n")
          break
        }

        # Update interval
        if(sign(diff_c$gamma_difference) == sign(diff_a$gamma_difference)) {
          a <- c
        } else {
          b <- c
        }

        iteration <- iteration + 1

        if(iteration == max_iterations) {
          cat("\nWarning: Maximum iterations (", max_iterations,
              ") reached without achieving desired tolerance.\n")
        }
      }

      # Calculate final point
      optimal_x <- (a + b) / 2
      optimal_y <- 2 * by_gP_mean - optimal_x

      gamma_result <- calculate_gamma_difference(by_data, by_gP_mean, optimal_x, gamma_actual)

      # Print final results for this By value
      cat("\nFinal Results for By =", by_val, ":\n")
      cat("Optimal point: (", sprintf("%.6f", optimal_x), ",", sprintf("%.6f", optimal_y), ")\n")
      cat("Final gamma achieved:", sprintf("%.2f", gamma_actual - gamma_result$gamma_difference), "%\n")
      cat("Points satisfying:", gamma_result$points_satisfying, "out of",
          gamma_result$total_points, "\n")
      cat("Number of iterations:", iteration, "\n")

      # Store results for this By value
      results[[as.character(by_val)]] <- list(
        By = by_val,
        g1 = optimal_x,
        g2 = optimal_y,
        gamma_difference = gamma_result$gamma_difference,
        gamma_achieved = gamma_actual - gamma_result$gamma_difference,
        points_satisfying = gamma_result$points_satisfying,
        total_points = gamma_result$total_points,
        iterations = iteration
      )
    }

    return(results)
  }

  update_ti_values <- function(data, optimal_results) {
    if (is.null(optimal_results)) return(data)

    for (by_val in names(optimal_results)) {
      result <- optimal_results[[by_val]]
      row_idx <- which(data$By == result$By)
      if (length(row_idx) > 0) {
        data$TI_Lower[row_idx] <- result$g1
        data$TI_Upper[row_idx] <- result$g2
      }
    }
    return(data)
  }

  ### plot cache
  interval_plots_cache <- reactiveVal(NULL)

  ### Observes
  observe({
    req(!reactivity_paused())
    if (!firstInputGiven() &&
        (length(input$fixed_var) > 0 ||
         length(input$random_var) > 0 ||
         input$multiplication_factor != 1 ||
         !is.null(input$by) ||
         input$tolerance_level != 50 ||
         input$confidence_of_tolerance != 95 ||
         input$percent_for_ci != 95 ||
         input$percent_for_pi != 95)) {
      firstInputGiven(TRUE)
      #updateTabsetPanel(session, "main_tabs", selected = "Formula")
    }
  }) #updates first-input-given condirion

  observe({
    req(!reactivity_paused())
    ConfidenceOfTolerance(input$confidence_of_tolerance)
  }) # updates conf of tolerance variable

  observe({
    req(input$multiplication_factor)
    req(!reactivity_paused())
    tryCatch({
      value <- as.numeric(input$multiplication_factor)
      if (is.na(value)) {
        stop("Multiplication Factor must be a numeric value")
      }
      MultiplicationFactorVariances(value)
    }, error = function(e) {
      ErrorMessage(as.character(e))
    })
  }) # updates error handling of multiplication factor

  observe({
    req(!reactivity_paused())
    req(input$fixed_var, input$random_var, input$multiplication_factor)
    # This will trigger the formulaString reactive to update
  }) # updates formula string rective

  observe({
    req(!reactivity_paused())
    input$fixed_var
    input$random_var
    input$multiplication_factor
    input$by
    input$tolerance_level
    input$confidence_of_tolerance
    input$percent_for_pi
    input$percent_for_ci
    input$covariate_1
    input$covariate_2
    input$covariate_3
    input$covariate_value_1
    input$covariate_value_2
    input$covariate_value_3

    if (intervalCalculationInitiated()) {
      input_changed(TRUE)
      last_change_time(Sys.time())
    }
  }) # monitors change in variables

  observe({
    req(!reactivity_paused())
    update_timer()
    if (intervalCalculationInitiated() && input_changed() && difftime(Sys.time(), last_change_time(), units = "secs") >= 1) {
      input_changed(FALSE)
      # Perform calculations without switching tabs
      withProgress(message = 'Updating calculations...', value = 0, {
        performCalculations(switch_tab = FALSE)
      })
    }
  }) # update if reactive timer is over a second

  observe({
    req(!reactivity_paused())
    req(IntervalsOutputInitial())

    # Scenario 1: IntervalsOutputInitial is not null but IntervalsOutputActual is null
    if (is.null(IntervalsOutputActual())) {
      lognorm_adjuster()
    }

    # Scenario 2: The input of the log_normal checkbox has changed
    observeEvent(input$log_normal, {
      lognorm_adjuster()
    }, ignoreInit = TRUE)

    # Scenario 3: The input of the show_reducible checkbox has changed
    observeEvent(input$show_reducible, {
      lognorm_adjuster()
    }, ignoreInit = TRUE)

    # Add new scenario for covariate changes
    observeEvent(input$covariate_var, {
      lognorm_adjuster()
    }, ignoreInit = TRUE)
  }) # observes lognorm input

  observe({
    req(!reactivity_paused())
    req(ForCDF1Dataframe())

    for (df in ForCDF1Dataframe()) {
      if (!is.null(df) && nrow(df) > 0) {
        local({
          local_df <- df
          by_value <- unique(local_df$By)

          # Initialize the plot (your existing plot code remains unchanged)
          output[[paste0("cdf_plot_", by_value)]] <- renderPlotly({
            req(input[[paste0("show_", by_value)]])

            tryCatch({
              x_value_input <- input[[paste0("x_value_", by_value)]]
              parsed_input <- parse_x_value_input(x_value_input)

              print(paste("Data for", by_value, ":"))
              print(str(local_df))
              print(paste("X value input:", x_value_input))

              # Check if this is a single chain dataset
              is_single_chain <- (single_chain_selected() && by_value == "1" || by_value == "SingleChain")

              # Set the title based on whether it's a single chain dataset
              plot_title <- if(is_single_chain) {
                "Posterior Predictive Distribution"
              } else {
                paste("Posterior Predictive Distribution for", by_value)
              }

              # Create the base plot using the display x values
              p <- plot_ly(data = local_df, x = ~x) %>%
                add_trace(y = ~y_mean, type = 'scatter', mode = 'lines',
                          name = "PP Distribution", line = list(color = 'rgba(0,100,80,1)', width = 2)) %>%
                add_trace(y = ~ci_lower, type = 'scatter', mode = 'lines',
                          name = "CI Lower", line = list(color = 'red', width = 2)) %>%
                add_trace(y = ~ci_upper, type = 'scatter', mode = 'lines',
                          name = "CI Upper", line = list(color = 'red', width = 2)) %>%
                add_ribbons(y = ~y_mean, ymin = ~ci_lower, ymax = ~ci_upper,
                            fillcolor = 'rgba(255,0,0,0.3)', line = list(color = 'rgba(255,0,0,0)'),
                            name = "CI Range", showlegend = FALSE)

              if (!is.null(parsed_input)) {
                if (parsed_input$type == "single") {
                  # Handle single value case (existing logic)
                  x_value <- parsed_input$value

                  # Convert input x value to log space if log_normal is TRUE
                  x_value_log <- if(input$log_normal) log(x_value) else x_value

                  # Check if the log value is in range
                  if (x_value_log >= local_df$xmin_log[1] && x_value_log <= local_df$xmax_log[1]) {
                    # Always use interpolation for precise values
                    y_value <- approx(local_df$x_log, local_df$y_mean, xout = x_value_log)$y
                    ci_lower <- approx(local_df$x_log, local_df$ci_lower, xout = x_value_log)$y
                    ci_upper <- approx(local_df$x_log, local_df$ci_upper, xout = x_value_log)$y

                    print("Interpolated values:")
                    print(data.frame(x = x_value, y = y_value, ci_lower = ci_lower, ci_upper = ci_upper))

                    # Add marker at the display x value
                    p <- p %>% add_trace(x = x_value, y = y_value,
                                         type = 'scatter', mode = 'markers',
                                         marker = list(color = 'black', size = 8),
                                         name = "Selected X Value",
                                         text = sprintf("x: %.4f<br>Y: %.4f<br>CI: [%.4f, %.4f]",
                                                        x_value, y_value, ci_lower, ci_upper),
                                         hoverinfo = 'text')

                    # Add horizontal dotted line
                    p <- p %>% add_segments(x = local_df$xmin[1], xend = x_value, y = y_value, yend = y_value,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)

                    # Add vertical dotted line
                    p <- p %>% add_segments(x = x_value, xend = x_value, y = 0, yend = y_value,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)

                  } else {
                    p <- p %>% add_annotations(
                      text = paste("x value", x_value, "is outside the valid range"),
                      x = 0.5, y = 1, xref = "paper", yref = "paper",
                      showarrow = FALSE, font = list(color = 'red')
                    )
                  }
                } else if (parsed_input$type == "range") {
                  # Handle range case
                  xmin <- parsed_input$xmin
                  xmax <- parsed_input$xmax

                  # Convert input x values to log space if log_normal is TRUE
                  xmin_log <- if(input$log_normal) log(xmin) else xmin
                  xmax_log <- if(input$log_normal) log(xmax) else xmax

                  # Check if both values are in range
                  if (xmin_log >= local_df$xmin_log[1] && xmin_log <= local_df$xmax_log[1] &&
                      xmax_log >= local_df$xmin_log[1] && xmax_log <= local_df$xmax_log[1]) {

                    # Calculate probabilities for both points
                    y_min <- approx(local_df$x_log, local_df$y_mean, xout = xmin_log)$y
                    y_max <- approx(local_df$x_log, local_df$y_mean, xout = xmax_log)$y

                    # Calculate CI values for both points
                    ci_lower_min <- approx(local_df$x_log, local_df$ci_lower, xout = xmin_log)$y
                    ci_upper_min <- approx(local_df$x_log, local_df$ci_upper, xout = xmin_log)$y
                    ci_lower_max <- approx(local_df$x_log, local_df$ci_lower, xout = xmax_log)$y
                    ci_upper_max <- approx(local_df$x_log, local_df$ci_upper, xout = xmax_log)$y

                    # Add markers for both points
                    p <- p %>% add_trace(x = c(xmin, xmax), y = c(y_min, y_max),
                                         type = 'scatter', mode = 'markers',
                                         marker = list(color = 'black', size = 8),
                                         name = "Range Bounds",
                                         text = c(
                                           sprintf("xmin: %.4f<br>P(X ≤ xmin): %.4f<br>CI: [%.4f, %.4f]",
                                                   xmin, y_min, ci_lower_min, ci_upper_min),
                                           sprintf("xmax: %.4f<br>P(X ≤ xmax): %.4f<br>CI: [%.4f, %.4f]",
                                                   xmax, y_max, ci_lower_max, ci_upper_max)
                                         ),
                                         hoverinfo = 'text')

                    # Add vertical lines for both points
                    p <- p %>% add_segments(x = xmin, xend = xmin, y = 0, yend = y_min,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)

                    p <- p %>% add_segments(x = xmax, xend = xmax, y = 0, yend = y_max,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)

                    # Add horizontal lines for both points
                    p <- p %>% add_segments(x = local_df$xmin[1], xend = xmin, y = y_min, yend = y_min,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)

                    p <- p %>% add_segments(x = local_df$xmin[1], xend = xmax, y = y_max, yend = y_max,
                                            line = list(color = 'black', width = 1, dash = 'dot'),
                                            showlegend = FALSE)
                  } else {
                    p <- p %>% add_annotations(
                      text = paste("One or both x values are outside the valid range"),
                      x = 0.5, y = 1, xref = "paper", yref = "paper",
                      showarrow = FALSE, font = list(color = 'red')
                    )
                  }
                }
              }

              p %>% layout(title = plot_title,
                           xaxis = list(title = "Value"),
                           yaxis = list(title = "Probability"))
            }, error = function(e) {
              print(paste("Error in plot for", by_value, ":", e$message))
              plot_ly() %>%
                add_annotations(
                  text = paste("Error in plot:", e$message),
                  showarrow = FALSE,
                  font = list(color = '#ff0000')
                )
            })
          })

          # Initialize the table with dynamic CI percentage in header
          output[[paste0("cdf_table_", by_value)]] <- renderDT({
            req(input[[paste0("show_", by_value)]])

            # Use isolate only for the CI percentage to prevent reactivity during rendering
            ci_percent <- isolate(input[[paste0("ci_", by_value)]])

            # Don't isolate x_value_input so it will trigger updates
            x_value_input <- input[[paste0("x_value_", by_value)]]

            # Parse the x value input with error handling
            parsed_input <- tryCatch({
              parse_x_value_input(x_value_input)
            }, error = function(e) {
              NULL
            })

            # Determine column names based on parsed input
            if (is.null(parsed_input)) {
              col_names <- c(
                "x value",
                "Probability of X ≤ x",
                paste0(ci_percent, "% CI")
              )
              is_range <- FALSE
            } else {
              is_range <- parsed_input$type == "range"
              if (is_range) {
                col_names <- c(
                  "x Interval",
                  "P(xmin ≤ X ≤ xmax)",
                  paste0(ci_percent, "% CI for P(xmin ≤ X ≤ xmax)")
                )
              } else {
                col_names <- c(
                  "x value",
                  "Probability of X ≤ x",
                  paste0(ci_percent, "% CI")
                )
              }
            }

            # Create the table data with fixed column names
            if (is.null(parsed_input)) {
              # Empty table case
              table_data <- data.frame(
                x_value = NA_character_,
                probability = NA_character_,
                ci_range = NA_character_,
                stringsAsFactors = FALSE
              )
            } else if (parsed_input$type == "single") {
              # Single value case
              x_value <- parsed_input$value
              x_value_log <- if(input$log_normal) log(x_value) else x_value

              if (x_value_log >= local_df$xmin_log[1] && x_value_log <= local_df$xmax_log[1]) {
                y_value <- approx(local_df$x_log, local_df$y_mean, xout = x_value_log)$y
                ci_lower <- approx(local_df$x_log, local_df$ci_lower, xout = x_value_log)$y
                ci_upper <- approx(local_df$x_log, local_df$ci_upper, xout = x_value_log)$y

                table_data <- data.frame(
                  x_value = as.character(round(x_value, 4)),
                  probability = as.character(round(y_value, 4)),
                  ci_range = sprintf("[%.4f, %.4f]", ci_lower, ci_upper),
                  stringsAsFactors = FALSE
                )
              } else {
                table_data <- data.frame(
                  x_value = as.character(x_value),
                  probability = "Out of range",
                  ci_range = "Out of range",
                  stringsAsFactors = FALSE
                )
              }
            } else if (parsed_input$type == "range") {
              # Range case
              xmin <- parsed_input$xmin
              xmax <- parsed_input$xmax

              xmin_log <- if(input$log_normal) log(xmin) else xmin
              xmax_log <- if(input$log_normal) log(xmax) else xmax

              if (xmin_log >= local_df$xmin_log[1] && xmin_log <= local_df$xmax_log[1] &&
                  xmax_log >= local_df$xmin_log[1] && xmax_log <= local_df$xmax_log[1]) {

                # Get the data for this by_value
                pred_dist <- predictedDistributions()
                pred_dist_by <- pred_dist[pred_dist$By == by_value, ]

                # Calculate CDF differences
                cdf_diffs <- sapply(1:nrow(pred_dist_by), function(i) {
                  mu <- pred_dist_by$median[i]
                  sigma <- sqrt(pred_dist_by$variance[i])
                  p_xmin <- pnorm(xmin_log, mean = mu, sd = sigma)
                  p_xmax <- pnorm(xmax_log, mean = mu, sd = sigma)
                  return(p_xmax - p_xmin)
                })

                # Calculate confidence interval
                alpha <- 1 - ci_percent/100
                ci_diff_lower <- quantile(cdf_diffs, alpha/2)
                ci_diff_upper <- quantile(cdf_diffs, 1 - alpha/2)

                # Calculate mean difference
                mean_diff <- mean(cdf_diffs)

                # Calculate display probabilities
                y_min <- approx(local_df$x_log, local_df$y_mean, xout = xmin_log)$y
                y_max <- approx(local_df$x_log, local_df$y_mean, xout = xmax_log)$y
                range_prob <- y_max - y_min

                table_data <- data.frame(
                  x_value = sprintf("[%.4f, %.4f]", xmin, xmax),
                  probability = as.character(round(range_prob, 4)),
                  ci_range = sprintf("[%.4f, %.4f]", ci_diff_lower, ci_diff_upper),
                  stringsAsFactors = FALSE
                )
              } else {
                table_data <- data.frame(
                  x_value = sprintf("[%.4f, %.4f]", xmin, xmax),
                  probability = "Out of range",
                  ci_range = "Out of range",
                  stringsAsFactors = FALSE
                )
              }
            }

            # Use the pre-determined column names for the datatable
            DT::datatable(
              table_data,
              options = list(
                pageLength = 10,
                dom = 't',
                ordering = FALSE,
                searching = FALSE,
                columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ),
              style = 'bootstrap4',
              rownames = FALSE,
              colnames = col_names
            )
          })

          # Add observer for CI slider
          observeEvent(input[[paste0("ci_", by_value)]], {
            if (table_rendering()) {
              # If table is rendering, revert the slider to its previous value
              updateSliderInput(session,
                                paste0("ci_", by_value),
                                value = isolate(input[[paste0("ci_", by_value)]]))
            }
          })
        })
      }
    }
  }) # dist plots generator

  observe({
    req(!reactivity_paused())

    # This will trigger whenever any of these inputs change
    input$show_reducible
    input$tolerance_level
    input$percent_for_pi

    # Check if show_reducible is TRUE and the probability levels don't match
    if (input$show_reducible && input$tolerance_level != input$percent_for_pi) {
      showNotification(
        "Warning: Tolerance Interval probability level and Prediction Interval percentage should be equal for valid reducible imprecision analysis.",
        type = "warning",
        duration = 10
      )
    }
  }) # error message regarding reducilbe uncertainty handler

  observe({
    req(!reactivity_paused())
    req(predictedDistributions())

    by_values <- unique(predictedDistributions()$By)

    for (by_value in by_values) {
      local({
        local_by_value <- by_value

        # Create a reactive observer for each x_value input
        observeEvent(input[[paste0("x_value_", local_by_value)]], {
          # Only trigger if the corresponding "show" checkbox is checked
          if (isTRUE(input[[paste0("show_", local_by_value)]])) {
            # Force the table to update by invalidating it
            table_id <- paste0("cdf_table_", local_by_value)

            # Set a flag to indicate table rendering is in progress
            table_rendering(TRUE)

            # Use invalidateLater with a very short delay to trigger re-rendering
            invalidateLater(10)

            # After a short delay, reset the rendering flag
            later::later(function() {
              table_rendering(FALSE)
            }, delay = 0.1)
          }
        }, ignoreInit = TRUE)
      })
    }
  }) # observer to update table when x value changes

  observe({
    # Check if any slider is at 100%
    if (input$tolerance_level == 100 ||
        input$confidence_of_tolerance == 100 ||
        input$percent_for_pi == 100 ||
        input$percent_for_ci == 100) {

      # Show a professional statistical warning
      showNotification(
        HTML("<b>Warning:</b> 100% probability intervals approach infinity in continuous distributions.
           For reliable results, please consider using 99% or 99.5% instead."),
        type = "warning",
        duration = 8
      )
    }
  }) # new warning interval sliders

  observe({
    # Get all By values from predictedDistributions
    req(predictedDistributions())
    by_values <- unique(predictedDistributions()$By)

    # Check each posterior predictive distribution CI slider
    for(by_value in by_values) {
      if (!is.null(input[[paste0("ci_", by_value)]]) &&
          input[[paste0("ci_", by_value)]] == 100) {

        showNotification(
          HTML("<b>Warning:</b> 100% credible intervals in the posterior predictive distribution
             cannot be computed exactly. The calculation will use a close approximation, but for more
             reliable visualization, consider using 99% or 99.9% instead."),
          type = "warning",
          duration = 10
        )
        break  # Show notification only once sliders are at 100%
      }
    }
  }) # new warning for posterior predictive CI sliders

  observe({
    # Reset optimal points when these inputs change
    input$fixed_var
    input$random_var
    input$multiplication_factor
    input$by
    input$tolerance_level
    input$confidence_of_tolerance
    input$covariate_1
    input$covariate_2
    input$covariate_3
    input$covariate_value_1
    input$covariate_value_2
    input$covariate_value_3

    # Reset the optimal points
    optimalPointResults(NULL)
  }) # Resets two sided TI optimal points when inputs change

  observe({
    req(firstTICalculation())
    req(!reactivity_paused())
    req(predictedDistributions())

    # Watch the sliders
    input$tolerance_level
    input$confidence_of_tolerance

    # Don't trigger on first run
    isolate({
      if (firstTICalculation()) {
        withProgress(message = 'Recalculating two-sided tolerance intervals...', value = 0, {
          # Calculate gP_mean
          scatter_data <- predictedDistributions()
          incProgress(0.3)

          # Calculate gP_mean and get the vector of means
          gP_mean_vector <- calculate_gP_mean(scatter_data, input$tolerance_level)
          incProgress(0.3)

          # Calculate optimal points using the confidence of tolerance level as gamma
          gamma_actual <- input$confidence_of_tolerance

          tryCatch({
            # Find optimal points for all By values
            results <- find_optimal_point(scatter_data, gP_mean_vector, gamma_actual)

            # Store both original and current values
            originalOptimalPointResults(results)

            # If log_normal is checked, transform the values
            if (input$log_normal) {
              transformed_results <- lapply(results, function(point) {
                if (!is.null(point$g1) && !is.null(point$g2)) {
                  result <- point
                  result$g1 <- exp(point$g1)
                  result$g2 <- exp(point$g2)
                  return(result)
                } else {
                  return(point)
                }
              })
              optimalPointResults(transformed_results)
            } else {
              optimalPointResults(results)
            }

          }, error = function(e) {
            print(paste("Error in find_optimal_point:", e$message))
            showNotification(paste("Error calculating optimal points:", e$message),
                             type = "error", duration = 10)
          })

          incProgress(0.4)
        })
      }
    })
  }) # Resets two sided TI optimal calculation when inputs change

  observe({
    # Reset optimal points when these inputs change
    input$fixed_var
    input$random_var
    input$multiplication_factor
    input$by
    input$tolerance_level
    input$confidence_of_tolerance
    input$covariate_var  # Add this line

    # Reset the optimal points
    optimalPointResults(NULL)
  }) # Checks status has krisnamoorthz TI been calculated yet

  observe({
    req(predictedDistributions())

    by_values <- unique(predictedDistributions()$By)
    states <- lapply(by_values, function(by_value) {
      input[[paste0("show_scatter_", by_value)]]
    })
    names(states) <- by_values

    scatter_checkbox_states(states)
  }) # records the scatter checkboxes in two sided intervals

  observe({
    req(input$krishnamoorthy_selected_by)
    krishnamoorthy_selected_by(input$krishnamoorthy_selected_by)
  }) # records the krishnamoorthy checkboxes in two sided intervals

  observe({
    req(InputData(), input$covariate_var)

    if(!is.null(input$covariate_var) && input$covariate_var != "") {
      # Check if the selected covariate column is numeric
      covariate_data <- InputData()[[input$covariate_var]]

      if(!is.numeric(covariate_data)) {
        showNotification(
          paste("Warning: Covariate column '", input$covariate_var,
                "' is not numeric. Please select a numeric column."),
          type = "warning",
          duration = 10
        )
      }

      # Check for NA values
      if(any(is.na(covariate_data))) {
        showNotification(
          paste("Warning: Covariate column '", input$covariate_var,
                "' contains NA values which may affect calculations."),
          type = "warning",
          duration = 10
        )
      }
    }
  }) # Covariate validation observer

  observe({
    req(InputData(), input$covariate_row)

    if(!is.null(input$covariate_row) && input$covariate_row != "None") {
      # Check if the selected row column is numeric
      row_data <- InputData()[[input$covariate_row]]
      if(!is.numeric(row_data)) {
        showNotification(
          paste("Warning: Covariate associated row '", input$covariate_row,
                "' is not numeric. Please select a numeric column."),
          type = "warning",
          duration = 10
        )
      }

      # Check for NA values
      if(any(is.na(row_data))) {
        showNotification(
          paste("Warning: Covariate associated row '", input$covariate_row,
                "' contains NA values which may affect calculations."),
          type = "warning",
          duration = 10
        )
      }
    }
  }) # Observer to validate covariate row selection


  ### Events
  observeEvent(input$load_data, {
    req(input$file)
    req(!reactivity_paused())
    # Initialize variables
    error_msg <- NULL
    IntervalsOutputInitial(NULL)
    IntervalsOutputActual(NULL)

    # Attempt to read the file
    data <- tryCatch({
      # First attempt: Try using fread (data.table's fast reader)
      data <- fread(input$file$datapath, showProgress = FALSE)

      if (is.null(data) || ncol(data) < 2) {
        # If fread fails, try read_csv
        data <- read_csv(input$file$datapath, show_col_types = FALSE)

        if (is.null(data) || ncol(data) < 2) {
          # If read_csv fails or returns insufficient data, try read.csv
          data <- read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE)

          if (is.null(data) || ncol(data) < 2) {
            # If read.csv fails, try read.table with more flexible options
            data <- read.table(input$file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

            if (is.null(data) || ncol(data) < 2) {
              # If read.table with tab separator fails, try with other common separators
              separators <- c(",", ";", "|", "/")
              for (sep in separators) {
                data <- tryCatch({
                  read.table(input$file$datapath, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
                }, error = function(e) NULL)
                if (!is.null(data) && ncol(data) >= 2) break
              }
            }
          }
        }
      }

      if (is.null(data) || ncol(data) < 2) {
        stop("The file could not be read or contains insufficient data.")
      }

      data
    }, error = function(e) {
      error_msg <<- paste("Error reading file:", e$message)
      NULL
    })

    # Check if data was successfully read
    if (is.null(data)) {
      ErrorMessage(error_msg)
      showNotification(error_msg, type = "error")
      return()
    }

    # Convert to data.table if it's not already
    if (!is.data.table(data)) {
      setDT(data)
    }

    # If the file is successfully read, update the InputData
    InputData(data)

    # Perform data inspection
    inspection_results <- tryCatch({
      inspect_data(data)
    }, error = function(e) {
      list(warning = paste("Error during data inspection:", e$message))
    })

    DataInspectionResults(inspection_results)

    # After successfully loading data, make sure you have ALL these lines:
    updateSelectInput(session, "fixed_var", choices = c("", names(data)))
    updateSelectInput(session, "random_var", choices = c("", names(data)))
    updateSelectInput(session, "covariate_var", choices = c("", names(data)))
    # Update the covariate row choices (initially disabled)
    # Update the three covariate inputs
    for(i in 1:3) {
      updateSelectInput(session, paste0("covariate_", i),
                        choices = c("None", names(data)),
                        selected = "None")
    }

    # Update the "by" choices, keeping "Single Chain Dataset" as the first option
    updateSelectInput(session, "by",
                      choices = c("Single Chain Dataset", names(data)),
                      selected = "Single Chain Dataset")

    # Switch to Data Preview tab
    updateTabsetPanel(session, "main_tabs", selected = "Data Preview")

    showNotification("File loaded successfully", type = "message")

    # If there are inspection warnings, show a notification
    if (length(unlist(inspection_results)) > 0) {
      showNotification("Data loaded with warnings. Check the Error tab.", type = "warning")
    }

    # Clear any previous error message
    ErrorMessage(NULL)
  })

  observeEvent(c(input$fixed_var, input$random_var, input$multiplication_factor, input$covariate_var), {
    req(input$fixed_var, input$random_var)
    req(!reactivity_paused())

    if (!FormulaGenerated()) {
      FormulaGenerated(TRUE)
      updateTabsetPanel(session, "main_tabs", selected = "Formula")
    }
  })

  observeEvent(input$calculate_ci, {
    req(!reactivity_paused())
    intervalCalculationInitiated(TRUE)
    withProgress(message = 'Calculating intervals...', value = 0, {
      performCalculations(switch_tab = TRUE)  # Switch tab when button is clicked
    })
  })

  observeEvent(input$plot_distribution, {
    req(!reactivity_paused())
    req(IntervalsOutputActual(), predictedDistributions())
    updateTabsetPanel(session, "main_tabs", selected = "Post. Pred. Distribution")
  })

  observeEvent(input$by, {
    # Set the flag to pause reactivity
    reactivity_paused(TRUE)

    # Store the current value of input$by
    current_by <- input$by
    single_chain_selected(identical(input$by, "Single Chain Dataset"))

    # Check if both "Single Chain Dataset" and other columns are selected
    if (length(current_by) > 1 && "Single Chain Dataset" %in% current_by) {
      # Show error notification
      showNotification("You cannot select both 'Single Chain Dataset' and other columns. Defaulting to 'Single Chain Dataset' only.",
                       type = "error", duration = 10)

      # Reset to just "Single Chain Dataset"
      updateSelectInput(session, "by",
                        choices = c("Single Chain Dataset", names(InputData())),
                        selected = "Single Chain Dataset")

      # Update help text
      byHelpTextValue("You can select multiple columns. They will be concatenated if more than one is selected.")
      originalByColumns(NULL)

      # Resume reactivity and exit early
      reactivity_paused(FALSE)
      return()
    }

    if (is.null(current_by) || length(current_by) == 0) {
      # If input$by is empty, set default
      byHelpTextValue("You can select multiple columns. They will be concatenated if more than one is selected.")
      originalByColumns(NULL)

      # Update to default
      updateSelectInput(session, "by",
                        choices = c("Single Chain Dataset", names(InputData())),
                        selected = "Single Chain Dataset")
    } else if (length(current_by) > 1) {
      # Store the original columns
      originalByColumns(current_by)

      # Create merged column directly here
      tryCatch({
        data <- InputData()

        # For data.table, use a different approach to extract multiple columns
        if (is.data.table(data)) {
          # Convert column names to symbols and use .SD to select them
          cols_to_merge <- lapply(current_by, as.name)
          merged_column <- data[, do.call(paste, c(.SD, sep = "_")), .SDcols = current_by]
        } else {
          # For data.frame, use the original approach
          merged_column <- apply(data[, current_by, drop = FALSE], 1, paste, collapse = "_")
        }

        # Create new column name
        new_column_name <- "Merged_By_Columns"

        # Add the new column to InputData
        data[[new_column_name]] <- merged_column
        InputData(data)

        # Print first 5 values of the new column for debugging
        print(paste("First 5 values of", new_column_name, ":"))
        print(head(merged_column, 5))

        # Update the select input with the new merged column
        updateSelectInput(session, "by",
                          choices = c("Single Chain Dataset", names(data)),
                          selected = new_column_name)

        # Update help text for multiple selections
        byHelpTextValue(paste("Concatenated By values consisting of", paste(current_by, collapse = ", "), "is selected."))
      }, error = function(e) {
        ErrorMessage(paste("Error merging columns:", e$message))
        showNotification(paste("Error merging columns:", e$message), type = "error")

        # Revert to default
        updateSelectInput(session, "by",
                          choices = c("Single Chain Dataset", names(InputData())),
                          selected = "Single Chain Dataset")
      })
    } else if (current_by == "Merged_By_Columns" && !is.null(originalByColumns())) {
      # If the merged column is selected, use the original columns for the help text
      byHelpTextValue(paste("Concatenated By values consisting of", paste(originalByColumns(), collapse = ", "), "is selected."))
    } else {
      # Single column selected, revert to original text
      byHelpTextValue("You can select multiple columns. They will be concatenated if more than one is selected.")
      originalByColumns(NULL)
    }

    # Update InputData if necessary
    if (identical(current_by, "Single Chain Dataset") && !"SingleChain" %in% names(InputData())) {
      data <- InputData()
      data$SingleChain <- 1
      InputData(data)
    }

    # Resume reactivity
    reactivity_paused(FALSE)
  }, ignoreInit = TRUE) # careful. pauses reactivity

  observeEvent(input$plot_quantiles, {
    req(!reactivity_paused())
    req(predictedDistributions())

    # Set the flag for first calculation and plot visibility
    firstTICalculation(TRUE)
    two_sided_ti_initiated(TRUE)

    withProgress(message = 'Calculating two-sided tolerance intervals...', value = 0, {
      # Calculate gP_mean
      scatter_data <- predictedDistributions()
      incProgress(0.3)

      # Calculate gP_mean and get the vector of means
      gP_mean_vector <- calculate_gP_mean(scatter_data, input$tolerance_level)
      incProgress(0.3)

      # Calculate optimal points using the confidence of tolerance level as gamma
      gamma_actual <- input$confidence_of_tolerance

      tryCatch({
        # Find optimal points for all By values
        results <- find_optimal_point(scatter_data, gP_mean_vector, gamma_actual)

        # Store both original and current values
        originalOptimalPointResults(results)

        # If log_normal is checked, transform the values
        if (input$log_normal) {
          transformed_results <- lapply(results, function(point) {
            if (!is.null(point$g1) && !is.null(point$g2)) {
              result <- point
              result$g1 <- exp(point$g1)
              result$g2 <- exp(point$g2)
              return(result)
            } else {
              return(point)
            }
          })
          optimalPointResults(transformed_results)
        } else {
          optimalPointResults(results)
        }

      }, error = function(e) {
        print(paste("Error in find_optimal_point:", e$message))
        showNotification(paste("Error calculating optimal points:", e$message),
                         type = "error", duration = 10)
      })

      incProgress(0.4)
    })

    # Switch to the Two-Sided Tolerance Interval tab
    updateTabsetPanel(session, "main_tabs", selected = "Two-Sided Intervals")
  }) # plot quantiles button observer

  observeEvent(input$by, {
    byIdentifierChanged(TRUE)
  })

  observeEvent(input$show_reducible, {
    req(!reactivity_paused())
    req(IntervalsOutputActual())

    print("Reducible imprecision checkbox changed - forcing plot update")

    # Force recalculation of plots
    interval_plots_cache(render_intervals_plot_one_sided(IntervalsOutputActual(), input))
  }, ignoreInit = TRUE)

  # Observer to validate covariate selections
  observeEvent(c(InputData(), input$covariate_1, input$covariate_2, input$covariate_3), {
    req(InputData())

    # Check all three covariate pairs
    for(i in 1:3) {
      covariate_col <- input[[paste0("covariate_", i)]]

      if(!is.null(covariate_col) && covariate_col != "None" && covariate_col %in% names(InputData())) {
        # Check if the selected covariate column is numeric
        covariate_data <- InputData()[[covariate_col]]

        if(!is.numeric(covariate_data)) {
          showNotification(
            paste("Warning: Covariate", i, "column '", covariate_col,
                  "' is not numeric. Please select a numeric column."),
            type = "warning",
            duration = 10
          )
        }

        # Check for NA values
        if(any(is.na(covariate_data))) {
          showNotification(
            paste("Warning: Covariate", i, "column '", covariate_col,
                  "' contains NA values which may affect calculations."),
            type = "warning",
            duration = 10
          )
        }
      }
    }
  }, ignoreInit = TRUE)

  # Observer to warn if covariate value is set but no covariate column selected
  observeEvent(c(input$covariate_value_1, input$covariate_value_2, input$covariate_value_3), {
    for(i in 1:3) {
      covariate_col <- input[[paste0("covariate_", i)]]
      covariate_val <- input[[paste0("covariate_value_", i)]]

      if(!is.null(covariate_val) && !is.na(covariate_val) && covariate_val != 1) {
        if(is.null(covariate_col) || covariate_col == "None") {
          showNotification(
            paste("Warning: Covariate Value", i, "is set to", covariate_val,
                  "but no Covariate", i, "column is selected. This value will be ignored."),
            type = "warning",
            duration = 5
          )
        }
      }
    }
  }, ignoreInit = TRUE)


  ### Outputs
  output$results_table <- renderDT({
    req(predictedDistributions())
    req(!reactivity_paused())
    print("Rendering PredictedDistributions")
    datatable(predictedDistributions(),
              options = list(pageLength = 10, dom = 't'),
              style = 'bootstrap4',
              rownames = FALSE
    )
  })

  output$ci_output <- renderDT({
    req(!reactivity_paused())
    req(IntervalsOutputActual())
    render_intervals_table(IntervalsOutputActual(), input, is_one_sided = TRUE)
  })

  output$intervals_plot_upper <- renderPlotly({
    req(!reactivity_paused())
    req(interval_plots_cache())

    print("Rendering upper intervals plot")
    print(paste("Show reducible:", input$show_reducible))  # Check checkbox state

    interval_plots_cache()$upper
  })

  output$intervals_plot_lower <- renderPlotly({
    req(!reactivity_paused())
    req(interval_plots_cache())

    print("Rendering lower intervals plot")
    print(paste("Show reducible:", input$show_reducible))  # Check checkbox state

    interval_plots_cache()$lower
  })

  output$quantile_ci_output <- renderDT({
    req(!reactivity_paused())
    req(IntervalsOutputActual())
    req(optimalPointResults())

    # Update TI values with Krishnamoorthy results
    updated_data <- update_ti_values(IntervalsOutputActual(), optimalPointResults())

    render_intervals_table(updated_data, input, is_one_sided = FALSE)
  })

  output$data_preview <- renderDT({
    req(!reactivity_paused())
    req(InputData())
    print("Rendering InputData")
    datatable(InputData(),
              options = list(
                pageLength = 10,
                dom = 'tp'  # 't' for table, 'p' for pagination
              ),
              style = 'bootstrap4',
              rownames = FALSE
    )
  })

  output$formula_output <- renderUI({
    req(!reactivity_paused())
    formulaString()
  })

  output$error_output <- renderText({
    error_msg <- ErrorMessage()
    if (!is.null(error_msg) && nchar(error_msg) > 0) {
      return(error_msg)
    } else {
      inspection_results <- DataInspectionResults()
      if (length(unlist(inspection_results)) > 0) {
        warnings <- unlist(inspection_results)
        return(paste("Possible data integrity issues:\n\n",
                     paste(warnings, collapse = "\n"),
                     "\n\nPlease review the data preview carefully.",
                     sep = ""))
      } else {
        return("No errors or warnings to display.")
      }
    }
  })

  output$distribution_plots <- renderUI({
    req(predictedDistributions())
    req(!reactivity_paused())

    # Control panel
    control_panel <- fluidRow(
      column(12,
             h4("Posterior Predictive Distribution Plot Controls"),
             uiOutput("by_value_controls")
      )
    )

    # Plot outputs
    plot_outputs <- uiOutput("cdf_plots")

    # Combine control panel and plots
    tagList(
      control_panel,
      plot_outputs
    )
  })

  output$by_value_controls <- renderUI({
    req(predictedDistributions())
    req(!reactivity_paused())

    by_values <- unique(predictedDistributions()$By)

    # Get the actual column name used for "by"
    by_column_name <- if(input$by == "Single Chain Dataset") {
      "By Identifier"
    } else if(input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
      paste(originalByColumns(), collapse = ", ")
    } else {
      input$by
    }

    lapply(by_values, function(by_value) {
      # Check if this is a single chain dataset
      is_single_chain <- (single_chain_selected() && by_value == "1" || by_value == "SingleChain")

      # Adjust the label based on whether it's a single chain dataset
      by_label <- if(is_single_chain) {
        "By Identifier:"
      } else {
        paste(by_column_name, ":", by_value)
      }

      fluidRow(
        column(width = 2,  # Column for "By" title
               strong(by_label)
        ),
        column(width = 2,  # Column for checkbox (wider for reliable rendering)
               div(class = "small-text-controls",
               checkboxInput(paste0("show_", by_value), "Show", value = FALSE, width = "100%")
               )
        ),
        column(width = 3,  # Column for slider
               div(class = "small-text-controls",
                   div(style = "display: flex; flex-direction: column; height: 100%;",
                       div(style = "flex-grow: 1;"),
                       sliderInput(paste0("ci_", by_value), "Credibility Interval (%)",
                                   min = 0, max = 100, value = 90, step = 1, width = "100%")
                   )
               )
        ),
        column(width = 3,  # Broader column for number of evaluation points
               div(class = "small-text-controls",
                   div(style = "display: flex; flex-direction: column; height: 100%;",
                       div(style = "flex-grow: 1;"),
                       numericInput(paste0("points_", by_value), "Number of Evaluation Points",
                                    value = 20, min = 5, max = 200, width = "100%")
                   )
               )
        ),
        column(width = 2,  # Column for X value - Changed from numericInput to textInput
               div(class = "small-text-controls",
                   div(style = "display: flex; flex-direction: column; height: 100%;",
                       div(style = "flex-grow: 1;"),
                       textInput(paste0("x_value_", by_value), "x", value = "", width = "100%"),
                       tags$small(style = "font-size: 0.7em; color: #aaa;",
                                  "Single value or [min,max]")
                   )
               )
        )
      )
    })
  })

  output$byHelpText <- renderUI({
    tags$p(
      style = "font-style: italic; font-size: 0.8em;",
      byHelpTextValue()
    )
  })

  output$cdf_plots <- renderUI({
    req(ForCDF1Dataframe())
    req(!reactivity_paused())

    plot_and_table_list <- lapply(ForCDF1Dataframe(), function(df) {
      if (!is.null(df) && nrow(df) > 0) {
        by_value <- unique(df$By)
        if (isTRUE(input[[paste0("show_", by_value)]])) {
          tagList(
            plotlyOutput(paste0("cdf_plot_", by_value), height = 400),
            DTOutput(paste0("cdf_table_", by_value)),
            hr()  # Add a horizontal line for separation
          )
        }
      }
    })

    do.call(tagList, plot_and_table_list)
  })

  output$quantile_scatter_controls <- renderUI({
    req(predictedDistributions())
    req(input$by)

    # Reset the change tracker
    byIdentifierChanged(FALSE)

    tryCatch({
      by_values <- unique(predictedDistributions()$By)

      # Get the actual column name used for "by"
      by_column_name <- if(input$by == "Single Chain Dataset") {
        "By Identifier"
      } else if(input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
        paste(originalByColumns(), collapse = ", ")
      } else {
        input$by
      }

      # Get current checkbox states
      current_states <- scatter_checkbox_states()

      div(
        class = "well",
        # Checkboxes section
        fluidRow(
          column(12,
                 div(
                   class = "checkbox-container",
                   style = "margin-bottom: 15px;", # Reduced margin
                   lapply(by_values, function(by_value) {
                     # Check if this is a single chain dataset
                     is_single_chain <- (single_chain_selected() && by_value == "1" || by_value == "SingleChain")

                     # Adjust the label based on whether it's a single chain dataset
                     by_label <- if(is_single_chain) {
                       "Single Chain Dataset"
                     } else {
                       paste(by_column_name, ":", by_value)
                     }

                     # Get previous state or default to FALSE
                     is_checked <- if (!is.null(current_states[[as.character(by_value)]])) {
                       current_states[[as.character(by_value)]]
                     } else {
                       FALSE
                     }

                     div(
                       class = "checkbox-item",
                       checkboxInput(
                         paste0("show_scatter_", by_value),
                         by_label,
                         value = is_checked  # Use preserved state
                       )
                     )
                   })
                 )
          )
        ),

        # Information section - directly after checkboxes with no extra containers
        div(
          style = "color: #ffffff; padding: 10px 0;",
          p("The tolerance interval confidence level slider adjusts the confidence level of the ellipses.
                Higher values create larger ellipses containing more data points."),
          tags$hr(style = "border-color: #668167; margin: 15px 0;"),
          p("Bivariate normal distribution is assumed for plotting purposes.
                Distribution of upper and lower tolerance interval may deviate from this assumption.")
        )
      )
    }, error = function(e) {
      # If there's an error, return a message
      div(
        p("Unable to generate controls. Please check your data and 'By' identifier selection."),
        p("Error: ", e$message)
      )
    })
  })

  output$krishnamoorthy_controls <- renderUI({
    req(predictedDistributions())
    req(input$by)

    tryCatch({
      by_values <- unique(predictedDistributions()$By)

      # Get the actual column name used for "by"
      by_column_name <- if(input$by == "Single Chain Dataset") {
        "By Identifier"
      } else if(input$by == "Merged_By_Columns" && !is.null(originalByColumns())) {
        paste(originalByColumns(), collapse = ", ")
      } else {
        input$by
      }

      # Create labels for the radio buttons
      radio_labels <- sapply(by_values, function(by_value) {
        if(single_chain_selected() && by_value == "1" || by_value == "SingleChain") {
          "Single Chain Dataset"
        } else {
          paste(by_column_name, ":", by_value)
        }
      })

      # Set names for the radio button choices
      radio_choices <- setNames(by_values, radio_labels)

      # Get the current selection
      current_selection <- krishnamoorthy_selected_by()

      div(
        class = "well",
        fluidRow(
          column(12,
                 div(
                   class = "checkbox-container",
                   style = "margin-bottom: 15px;",
                   # Enhanced CSS styling
                   tags$style(HTML("
                            .krishnamoorthy-radio .radio-inline {
                                display: inline-block;
                                margin: 5px 10px;
                                background-color: #000000;
                                padding: 8px 15px;
                                border-radius: 4px;
                                border: 1px solid #668167;
                                min-width: 200px;
                                position: relative;
                            }
                            .krishnamoorthy-radio .radio-inline input[type='radio'] {
                                -webkit-appearance: none;
                                -moz-appearance: none;
                                appearance: none;
                                width: 14px;
                                height: 14px;
                                border: 2px solid #668167;
                                border-radius: 4px;
                                outline: none;
                                transition: background-color 0.3s;
                                vertical-align: middle;
                                margin-right: 8px;
                                margin-top: 0;
                                position: relative;
                                top: 1px;
                            }
                            .krishnamoorthy-radio .radio-inline input[type='radio']:checked {
                                background-color: #668167;
                            }
                            .krishnamoorthy-radio .radio-inline input[type='radio']:checked::after {
                                content: '';
                                display: block;
                                width: 3px;
                                height: 7px;
                                border: solid white;
                                border-width: 0 2px 2px 0;
                                transform: rotate(45deg);
                                margin: 1px 0 0 4px;
                            }
                            .krishnamoorthy-radio .radio-inline span {
                                vertical-align: middle;
                                margin-left: 4px;
                            }
                            /* Flex container for better alignment */
                            .krishnamoorthy-radio .radio-group {
                                display: flex;
                                flex-wrap: wrap;
                                gap: 10px;
                                justify-content: flex-start;
                            }
                        ")),
                   div(
                     class = "krishnamoorthy-radio",
                     div(
                       class = "radio-group",
                       radioButtons(
                         "krishnamoorthy_selected_by",
                         label = NULL,
                         choices = radio_choices,
                         selected = current_selection,  # Use the stored selection
                         inline = TRUE
                       )
                     )
                   )
                 )
          )
        )
      )
    }, error = function(e) {
      div(
        p("Unable to generate controls. Please check your data and 'By' identifier selection."),
        p("Error: ", e$message)
      )
    })
  })

  output$quantile_scatter_plot <- renderPlotly({
    req(two_sided_ti_initiated())
    req(predictedDistributions())
    req(input$by)
    req(!byIdentifierChanged())

    # Debounce the ellipse confidence input
    ellipse_level <- reactive({
      input$confidence_of_tolerance / 100
    }) %>% debounce(1000)

    tryCatch({
      # Get the data and ensure it's a data frame
      scatter_data <- as.data.frame(predictedDistributions())

      # Verify required columns exist
      req(all(c("lower_quantile", "upper_quantile", "By") %in% names(scatter_data)))

      # Calculate the actual percentiles for labels
      tolerance_level <- input$tolerance_level
      lower_percentile <- (100 - tolerance_level) / 2
      upper_percentile <- 100 - lower_percentile

      # Check if log_normal is checked
      is_log_normal <- input$log_normal

      # Initialize the plot
      p <- plot_ly()

      # Get unique By values
      by_values <- unique(scatter_data$By)

      # Create a color palette for the different By values
      n_colors <- length(by_values)
      colors <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"))(n_colors)

      # Performance improvement: Pre-calculate all ellipses before adding to plot
      ellipse_data <- list()

      # First pass: calculate all ellipses and means for checked datasets
      for(i in seq_along(by_values)) {
        by_value <- by_values[i]
        if(isTRUE(input[[paste0("show_scatter_", by_value)]])) {
          data_subset <- scatter_data[scatter_data$By == by_value, ]

          # Calculate mean and covariance for the current subset
          subset_data <- data_subset[, c("lower_quantile", "upper_quantile")]
          data_mean <- colMeans(subset_data)
          data_cov <- cov(subset_data)

          # Generate ellipse points with the user-specified confidence level
          ellipse_points <- create_confidence_ellipse_points(
            cov_matrix = data_cov,
            centre = data_mean,
            level = ellipse_level(),
            npoints = 100
          )

          # Store for later use
          ellipse_data[[by_value]] <- list(
            points = ellipse_points,
            mean = data_mean,
            subset = data_subset,
            is_single_chain = (single_chain_selected() && by_value == "1" || by_value == "SingleChain")
          )
        }
      }

      # Second pass: Add scatter points and mean points
      for(i in seq_along(by_values)) {
        by_value <- by_values[i]
        if(isTRUE(input[[paste0("show_scatter_", by_value)]]) && !is.null(ellipse_data[[by_value]])) {
          # Get stored data
          e_data <- ellipse_data[[by_value]]

          # Set the trace name
          trace_name <- if(e_data$is_single_chain) {
            "Single Chain Dataset"
          } else {
            paste("By:", by_value)
          }

          # Add scatter points (conditionally exponentiated)
          p <- add_trace(p,
                         data = e_data$subset,
                         x = if(is_log_normal) ~exp(upper_quantile) else ~upper_quantile,
                         y = if(is_log_normal) ~exp(lower_quantile) else ~lower_quantile,
                         type = 'scatter',
                         mode = 'markers',
                         name = trace_name,
                         marker = list(
                           size = 8,
                           color = colors[i],
                           opacity = 0.6
                         ),
                         hoverinfo = 'text',
                         text = ~paste(
                           "By:", by_value,
                           "<br>Lower Quantile:", round(if(is_log_normal) exp(lower_quantile) else lower_quantile, 4),
                           "<br>Upper Quantile:", round(if(is_log_normal) exp(upper_quantile) else upper_quantile, 4)
                         ),
                         showlegend = TRUE
          )

          # Add mean point (conditionally exponentiated)
          p <- add_trace(p,
                         x = if(is_log_normal) exp(e_data$mean[2]) else e_data$mean[2],
                         y = if(is_log_normal) exp(e_data$mean[1]) else e_data$mean[1],
                         type = 'scatter',
                         mode = 'markers',
                         name = paste(trace_name, "Mean"),
                         marker = list(
                           size = 12,
                           color = colors[i],
                           symbol = 'diamond'
                         ),
                         hoverinfo = 'text',
                         text = paste(
                           "Mean for", trace_name,
                           "<br>Lower Quantile:", round(if(is_log_normal) exp(e_data$mean[1]) else e_data$mean[1], 4),
                           "<br>Upper Quantile:", round(if(is_log_normal) exp(e_data$mean[2]) else e_data$mean[2], 4)
                         ),
                         showlegend = FALSE
          )
        }
      }

      # Third pass: Add all ellipses last (to keep them at the forefront)
      for(i in seq_along(by_values)) {
        by_value <- by_values[i]
        if(isTRUE(input[[paste0("show_scatter_", by_value)]]) && !is.null(ellipse_data[[by_value]])) {
          # Get stored data
          e_data <- ellipse_data[[by_value]]

          # Set the trace name
          trace_name <- if(e_data$is_single_chain) {
            "Single Chain Dataset"
          } else {
            paste("By:", by_value)
          }

          # Add ellipse with full opacity (conditionally exponentiated)
          base_color <- colors[i]
          darker_color <- colorRampPalette(c(base_color, "black"))(3)[2] # Mix with black for darker shade

          p <- add_trace(p,
                         x = if(is_log_normal) exp(e_data$points[, "y"]) else e_data$points[, "y"],
                         y = if(is_log_normal) exp(e_data$points[, "x"]) else e_data$points[, "x"],
                         type = 'scatter',
                         mode = 'lines',
                         name = paste(trace_name, "Ellipse"),
                         line = list(
                           color = darker_color,
                           width = 2.5
                         ),
                         hoverinfo = 'text',
                         text = paste(input$confidence_of_tolerance, "% Confidence Ellipse for", trace_name),
                         showlegend = TRUE
          )
        }
      }

      # Set layout with updated titles and labels (no [Exponentiated] text)
      p %>% layout(
        title = sprintf("Two-Sided Tolerance Intervals (%.1f%% Probability Level, %.1f%% Confidence Level)",
                        tolerance_level, input$confidence_of_tolerance),
        xaxis = list(
          title = sprintf("Upper Tolerance Interval Probability Level (%.1f%%)", upper_percentile),
          zeroline = TRUE
        ),
        yaxis = list(
          title = sprintf("Lower Tolerance Interval Probability Level (%.1f%%)", lower_percentile),
          zeroline = TRUE
        ),
        hovermode = 'closest',
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          y = -0.2
        )
      )
    }, error = function(e) {
      # If there's an error, return an empty plot with an error message
      plot_ly() %>%
        add_annotations(
          text = paste("Error in generating plot:", e$message,
                       "<br>Please check your data and 'By' identifier selection."),
          showarrow = FALSE,
          font = list(color = 'red')
        )
    })
  })

  output$krishnamoorthy_plot <- renderPlotly({
    req(two_sided_ti_initiated())
    req(predictedDistributions())
    req(input$by)
    req(!byIdentifierChanged())
    req(optimalPointResults())
    req(input$tolerance_level)
    req(input$confidence_of_tolerance)

    tryCatch({
      scatter_data <- as.data.frame(predictedDistributions())
      req(all(c("lower_quantile", "upper_quantile", "By") %in% names(scatter_data)))

      tolerance_level <- input$tolerance_level
      lower_percentile <- (100 - tolerance_level) / 2
      upper_percentile <- 100 - lower_percentile

      p <- plot_ly()

      selected_by <- krishnamoorthy_selected_by()

      if(is.null(selected_by) || selected_by == "") {
        return(plot_ly() %>% layout(title = "Two-Sided Tolerance Intervals"))
      }

      # Get the optimal points for the selected By value
      optimal_result <- optimalPointResults()[[as.character(selected_by)]]

      # Get the original data subset
      data_subset <- scatter_data[scatter_data$By == selected_by, ]
      is_single_chain <- single_chain_selected() && (selected_by == "1" || selected_by == "SingleChain")
      trace_name <- if(is_single_chain) "Single Chain Dataset" else paste("By:", selected_by)

      # Determine if we're in log-normal mode
      is_log_normal <- input$log_normal

      # Add scatter points
      p <- add_trace(p,
                     data = data_subset,
                     x = if(is_log_normal) ~exp(upper_quantile) else ~upper_quantile,
                     y = if(is_log_normal) ~exp(lower_quantile) else ~lower_quantile,
                     type = 'scatter',
                     mode = 'markers',
                     name = trace_name,
                     marker = list(
                       size = 8,
                       color = '#1f77b4',
                       opacity = 0.8
                     ),
                     hoverinfo = 'text',
                     text = ~paste(
                       "By:", selected_by,
                       "<br>Lower Quantile:", round(if(is_log_normal) exp(lower_quantile) else lower_quantile, 4),
                       "<br>Upper Quantile:", round(if(is_log_normal) exp(upper_quantile) else upper_quantile, 4)
                     ),
                     showlegend = TRUE
      )

      if(!is.null(optimal_result)) {
        # Get the original log-space values
        original_optimal_result <- isolate(originalOptimalPointResults())[[as.character(selected_by)]]
        g1_log <- original_optimal_result$g1
        g2_log <- original_optimal_result$g2
        gP_mean_log <- (g1_log + g2_log) / 2

        # Determine display values based on log_normal setting
        g1_display <- if(is_log_normal) exp(g1_log) else g1_log
        g2_display <- if(is_log_normal) exp(g2_log) else g2_log

        # Calculate the line
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

        if(is_log_normal) {
          # For log-normal case: y_new = exp(a*log(x_new)+b) where a=-1 and b=2*gP_mean
          x_line_display <- exp(x_line_log)
          y_line_display <- exp(-log(x_line_display) + 2*gP_mean_log)
        } else {
          # For non-log-normal case: y = 2*gP_mean - x in log space
          x_line_display <- x_line_log
          y_line_display <- 2 * gP_mean_log - x_line_log
        }

        # Add the Krishnamoorthy line
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
                       text = paste("Krishnamoorthy Line")
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
                         "Two-sided TI (Krishnamoorthy)<br>Lower TI Limit: %.4f<br>Upper TI Limit: %.4f<br>Achieved Value of Bisection Method: %.2f%%",
                         g1_display, g2_display, optimal_result$gamma_achieved
                       )
        )

        # Calculate display ranges for shapes
        x_min_display <- if(is_log_normal) exp(min(x_range)) else min(x_range)
        x_max_display <- if(is_log_normal) exp(max(x_range)) else max(x_range)
        y_min_display <- if(is_log_normal) exp(min(y_range)) else min(y_range)
        y_max_display <- if(is_log_normal) exp(max(y_range)) else max(y_range)

        # Calculate padding for display ranges
        x_padding <- (x_max_display - x_min_display) * 0.1
        y_padding <- (y_max_display - y_min_display) * 0.1

        # Set layout with shapes for lines and shaded area
        p <- p %>% layout(
          title = sprintf("Two-Sided Tolerance Intervals (%.1f%% Probability Level, %.1f%% Confidence Level)",
                          tolerance_level, input$confidence_of_tolerance),
          xaxis = list(
            title = sprintf("Upper Tolerance Interval Probability Level (%.1f%%)%s",
                            upper_percentile, if(is_log_normal) " [Exponentiated]" else ""),
            zeroline = TRUE
          ),
          yaxis = list(
            title = sprintf("Lower Tolerance Interval Probability Level (%.1f%%)%s",
                            lower_percentile, if(is_log_normal) " [Exponentiated]" else ""),
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

      # If no optimal results, return basic plot
      p %>% layout(
        title = sprintf("Two-Sided Tolerance Intervals (%.1f%% Probability Level, %.1f%% Confidence Level)",
                        tolerance_level, input$confidence_of_tolerance),
        xaxis = list(
          title = sprintf("Upper Tolerance Interval Probability Level (%.1f%%)%s",
                          upper_percentile, if(is_log_normal) " [Exponentiated]" else ""),
          zeroline = TRUE
        ),
        yaxis = list(
          title = sprintf("Lower Tolerance Interval Probability Level (%.1f%%)%s",
                          lower_percentile, if(is_log_normal) " [Exponentiated]" else ""),
          zeroline = TRUE
        ),
        hovermode = 'closest',
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          y = -0.2
        )
      )
    }, error = function(e) {
      plot_ly() %>%
        add_annotations(
          text = paste("Error in generating plot:", e$message,
                       "<br>Please check your data and 'By' identifier selection."),
          showarrow = FALSE,
          font = list(color = 'red')
        )
    })
  })

}



shinyApp(ui, server)
