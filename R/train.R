#'
#' Train a KNN-1 TWDTW model with optional GAM resampling
#'
#' This function prepares a KNN-1 model with the Time Warp Dynamic Time Warping (TWDTW) algorithm.
#' If a formula is provided, the training samples are resampled using Generalized Additive Models (GAM).
#'
#' @param x A two-dimensional stars object (x, y) with time and bands as attributes.
#' @param y An sf object with the coordinates of the training points.
#' @param formula Either NULL or a formula to reduce samples of the same label using Generalized Additive Models (GAM).
#' Default is \code{band ~ s(time)}. See details.
#' @param start_column Name of the column in y that indicates the start date. Default is 'start_date'.
#' @param end_column Name of the column in y that indicates the end date. Default is 'end_date'.
#' @param label_colum Name of the column in y containing land use labels. Default is 'label'.
#' @param sampling_freq The time frequency for sampling, including the unit (e.g., '16 day').
#' If NULL, the function will infer the frequency. This parameter is only used if a formula is provided.
#' @param ... Additional arguments passed to the GAM function.
#'
#' @details If \code{formula} is NULL, the KNN-1 model will retain all training samples. If a formula is passed (e.g., \code{band ~ \link[mgcv]{s}(time)}),
#' then samples of the same label (land cover class) will be resampled using GAM.
#' Resampling can significantly reduce prediction processing time.
#'
#' @return A 'knn1_twdtw' model containing the trained model information and the data used.
#'
#' @examples
#'
#' # Example will be added once the function is properly implemented and tested.
#' 
#' @export
knn1_twdtw <- function(x, y, formula = NULL, start_column = 'start_date',
                       end_column = 'end_date', label_colum = 'label',
                       sampling_freq = NULL, ...){

  # Check if x is a stars object with a time dimension
  if (!inherits(x, "stars") || length(dim(x)) != 2) {
    stop("x must be a stars object with two dimensions")
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
  y[, start_column] <- to_date_time(y[[start_column]])
  y[, end_column] <- to_date_time(y[[end_column]])

  # Prepare time series samples from stars object
  ts_data <- st_extract(x, y)
  ts_data$label <- y[[label_colum]]
  ts_data <- prepare_time_series(as.data.frame(ts_data))
  ts_data$ts_id <- NULL

  if(!is.null(formula)) {

    # Check if formula has two
    if(length(all.vars(formula)) != 2) {
      stop("The formula should have only one predictor!")
    }

    # Determine sampling frequency
    if (is.null(sampling_freq)) {
      sampling_freq <- get_time_series_freq(ts_data)
    }

    # Shift dates
    ts_data$observations <- lapply(ts_data$observations, shift_ts_dates)

    # Split data frame by label
    ts_data <- unnest(ts_data, cols = .data$observations)
    ts_data <- nest(ts_data, .by = .data$label, .key = "observations")

    # Define GAM function
    gam_fun <- function(band, t, pred_t, formula, ...){
      df <- setNames(list(band, as.numeric(t)), all.vars(formula))
      pred_t[[all.vars(formula)[2]]] <- as.numeric(pred_t[[all.vars(formula)[2]]])
      fit <- mgcv::gam(data = df, formula = formula, ...)
      predict(fit, newdata = pred_t)
    }

    # Apply GAM function
    ts_data$observations <- lapply(ts_data$observations, function(ts){
      y_time <- ts$time
      ts$time <- NULL
      pred_time <- setNames(list(seq(min(y_time), max(y_time), by = sampling_freq)), all.vars(formula)[2])
      cbind(pred_time, as.data.frame(sapply(as.list(ts), function(band) {
        gam_fun(band, y_time, pred_time, formula, ...)
      })))
    })

  }

  model <- list()
  model$call <- match.call()
  model$formula <- formula
  model$data <- ts_data
  class(model) <- "knn1_twdtw"

  return(model)

}

#' Compute the Most Common Sampling Frequency in a Stars Object
#'
#' This function calculates the most common difference between consecutive time points in a stars object. 
#' This can be useful for determining the sampling frequency of the time series data.
#'
#' @param x A stars object containing time series data.
#'
#' @return A difftime object representing the most common time difference between consecutive samples.
#'
#'
get_stars_time_freq <- function(x) {

  # Extract the time dimension
  time_values <- st_get_dimension_values(x, "time")

  # Compute the differences between consecutive time points
  time_diffs <- diff(time_values)

  # Convert differences to days (while retaining the difftime class)
  time_diffs <- as.difftime(time_diffs, units = "days")

  # Identify the mode
  mode_val_index <- which.max(tabulate(match(time_diffs, unique(time_diffs))))
  freq <- diff(time_values[mode_val_index:(mode_val_index+1)])

  return(freq)
}

get_time_series_freq <- function(x) {

  # Extract the time dimension
  time_values <- unlist(lapply(x$observations, function(ts) ts$time))

  # Compute the differences between consecutive time points
  time_diffs <- unlist(lapply(x$observations, function(ts) diff(ts$time)))

  # Convert differences to days (while retaining the difftime class)
  time_diffs <- as.difftime(time_diffs, units = "days")

  # Identify the mode
  mode_val_index <- which.max(tabulate(match(time_diffs, unique(time_diffs))))
  freq <- diff(time_values[mode_val_index:(mode_val_index+1)])

  return(freq)

}