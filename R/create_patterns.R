#' Create a Pattern Using GAM
#'
#' This function creates a pattern based on Generalized Additive Models (GAM).
#' It uses the specified formula to fit the model and predict values.
#' It operates on stars objects (for spatial-temporal data) and sf objects (for spatial data).
#'
#' @param x A stars object representing spatial-temporal data.
#' @param y An sf object representing spatial data with associated attributes.
#' @param formula A formula for the GAM. 
#' @param start_column Name of the column in y that indicates the start date. Default is 'start_date'.
#' @param end_column Name of the column in y that indicates the end date. Default is 'end_date'.
#' @param label_colum Name of the column in y that contains labels. Default is 'label'.
#' @param sampling_freq The time frequency for sampling. If NULL, the function will infer it.
#' @param ... Additional arguments passed to the GAM function.
#'
#' @return A list containing the predicted values for each label.
#'
#' 
#'
#' @export
create_patterns = function(x, y, formula, start_column = 'start_date', 
                          end_column = 'end_date', label_colum = 'label', 
                          sampling_freq = NULL, ...){
  
  # Check if x is a stars object with a time dimension
  if (!inherits(x, "stars") || dim(x)['time'] < 1) {
    stop("x must be a stars object with a 'time' dimension")
  }

  # Check if y is an sf object with point geometry
  if (!inherits(y, "sf") || !all(st_is(y, "POINT"))) {
    stop("y must be an sf object with point geometry")
  }

  # Check for required columns in y
  required_columns <- c(start_column, end_column, label_colum)
  missing_columns <- setdiff(required_columns, names(y))
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns in y:", paste(missing_columns, collapse = ", ")))
  }

  # Convert columns to date-time
  y[ , start_column] <- to_date_time(y[[start_column]])
  y[ , end_column] <- to_date_time(y[[end_column]])

  # Extract time series from stars
  y_ts <- extract_time_series(x, y)
  y_ts$geom <- NULL

  # Shift dates
  unique_ids <- unique(y_ts$id)
  for (id in seq_along(unique_ids)) {
    idx <- y_ts$id == unique_ids[id]
    y_ts[idx, 'time'] <- shift_dates(y_ts[idx, 'time'])
  }

  # Split data frame by label and remove label column
  y_ts <- lapply(split(y_ts, y_ts$label), function(df) {
    df$label <- NULL
    return(df)
  })

  # Determine sampling frequency
  if (is.null(sampling_freq)) {
    warning("Sampling frequency not defined; inferring from the stars object.")
    sampling_freq <- get_stars_time_freq(x)
  }

  # Extract variables from formula
  vars <- all.vars(formula)

  # Define GAM function
  gam_fun = function(y, t, formula_vars, ...){
    df = data.frame(y, t = as.numeric(t))
    names(df) <- formula_vars
    fit = mgcv::gam(data = df, formula = as.formula(paste(formula_vars[1], "~", formula_vars[2])), ...)
    pred_t = data.frame(t = as.numeric(seq(min(t), max(t), by = sampling_freq)))
    predict(fit, newdata = pred_t)
  }

  # Apply GAM function
  res <- lapply(y_ts, function(ts){
    y_time <- ts$time
    ts$time <- NULL
    ts$id <- NULL
    sapply(as.list(ts), function(y) {
      gam_fun(y, y_time, vars)
    })
  })

  return(res)
}



#' Extract Time Series from a Stars Object for Specified Points
#'
#' This function extracts a time series from a stars object for each specified point in the sf object. 
#' Each extracted sample is then labeled with an ID and the label from the sf object.
#'
#' @param x A stars object containing the raster time series data.
#' @param y An sf object containing the point geometries and their associated labels.
#'
#' @return A data.frame with the extracted time series for each point in the sf object, 
#'         with additional columns for the ID and label of each sample.
#'
#'
extract_time_series <- function(x, y) {
  ts_samples <- st_extract(x, y)
  ts_samples$id <- 1:dim(ts_samples)["geom"]
  ts_samples$label <- y$label
  as.data.frame(ts_samples)
}



#' Compute the Most Common Sampling Frequency in a Stars Object
#'
#' This function calculates the most common difference between consecutive time points in a stars object. 
#' This can be useful for determining the sampling frequency of the time series data.
#'
#' @param stars_object A stars object containing time series data.
#'
#' @return A difftime object representing the most common time difference between consecutive samples.
#'
#'
get_stars_time_freq <- function(stars_object) {

  # Extract the time dimension
  time_values <- st_get_dimension_values(stars_object, "time")

  # Compute the differences between consecutive time points
  time_diffs <- diff(time_values)

  # Convert differences to days
  time_diffs <- as.difftime(time_diffs, units = "days")

  # Compute the most common difference (mode)
  freq <- unique(time_diffs)[which.max(tabulate(match(time_diffs, unique(time_diffs))))]

  return(freq)
}
