#'
#' Train a KNN-1 TWDTW model with optional GAM resampling
#'
#' This function prepares a KNN-1 model with the Time Warp Dynamic Time Warping (TWDTW) algorithm.
#' If a formula is provided, the training samples are resampled using Generalized Additive Models (GAM).
#'
#' @param x A three-dimensional stars object (x, y, time) with bands as attributes.
#' @param y An sf object with the coordinates of the training points.
#' @param time_weight A numeric vector with length two (steepness and midpoint of logistic weight) or a function.
#' See details in \link[twdtw]{twdtw}.
#' @param cycle_length The length of the cycle, e.g. phenological cycles. Details in \link[twdtw]{twdtw}.
#' @param time_scale Specifies the time scale for the observations. Details in \link[twdtw]{twdtw}.
#' @param formula Either NULL or a formula to reduce samples of the same label using Generalized Additive Models (GAM).
#' Default is \code{band ~ s(time)}. See details.
#' @param start_column Name of the column in y that indicates the start date. Default is 'start_date'.
#' @param end_column Name of the column in y that indicates the end date. Default is 'end_date'.
#' @param label_colum Name of the column in y containing land use labels. Default is 'label'.
#' @param sampling_freq The time frequency for sampling, including the unit (e.g., '16 day').
#' If NULL, the function will infer the frequency. This parameter is only used if a formula is provided.
#' @param ... Additional arguments passed to the \link[mgcv]{gam} function and to \link[twdtw]{twdtw} function.
#'
#' @details If \code{formula} is NULL, the KNN-1 model will retain all training samples. If a formula is passed (e.g., \code{band ~ \link[mgcv]{s}(time)}),
#' then samples of the same label (land cover class) will be resampled using GAM.
#' Resampling can significantly reduce prediction processing time.
#'
#' @return A 'twdtw_knn1' model containing the trained model information and the data used.
#'
#' @examples
#' \dontrun{
#'
#' # Read training samples
#' samples_path <- 
#'   system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat")
#' 
#' samples <- st_read(samples_path, quiet = TRUE)
#'
#' # Get satellite image time sereis files
#' tif_path <- system.file("mato_grosso_brazil", package = "dtwSat")
#' tif_files <- dir(tif_path, pattern = "\\.tif$", full.names = TRUE)
#'
#' # Get acquisition dates
#' acquisition_date <- regmatches(tif_files, regexpr("[0-9]{8}", tif_files))
#' acquisition_date <- as.Date(acquisition_date, format = "%Y%m%d")
#'
#' # Create a 3D datacube
#' dc <- read_stars(tif_files,
#'                  proxy = FALSE,
#'                  along = list(time = acquisition_date),
#'                  RasterIO = list(bands = 1:6))
#' dc <- st_set_dimensions(dc, 3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR"))
#' dc <- split(dc, c("band"))
#'
#' # Create a knn1-twdtw model
#' m <- twdtw_knn1(x = dc,
#'                 y = samples,
#'                 cycle_length = 'year',
#'                 time_scale = 'day',
#'                 time_weight = c(steepness = 0.1, midpoint = 50),
#'                 formula = band ~ s(time))
#'
#' print(m)
#'
#' # Visualize model patterns
#' plot(m)
#'
#' # Classify satellite images
#' system.time(lu <- predict(dc, model = m))
#'
#' # Visualise land use classification
#' ggplot() +
#'   geom_stars(data = lu) +
#'   theme_minimal()
#'
#' }
#' @export
twdtw_knn1 <- function(x, y, time_weight, cycle_length, time_scale,
                       formula = NULL, start_column = 'start_date',
                       end_column = 'end_date', label_colum = 'label',
                       sampling_freq = NULL, ...){

  # Check if x is a stars object with a time dimension
  if (!inherits(x, "stars") || dim(x)['time'] < 1 || length(dim(x)) != 3) {
    stop("x must be a three-dimensional stars object with a 'time' dimension")
  }

  x <- split(x, c("time"))

  # Check if y is an sf object with point geometry
  if (!inherits(y, "sf") || !all(st_is(y, "POINT"))) {
    stop("y must be an sf object with point geometry")
  }

  # check for minimum set of twdtw arguments
  if (!(is.function(time_weight) || (is.numeric(time_weight) && length(time_weight) == 2))) stop("'time_weight' should be either a function or a numeric vector with length two")
  if (is.null(cycle_length)) stop("The 'cycle_length' argument is missing.")
  if (is.null(time_scale)) stop("The 'time_scale' argument is missing for 'cycle_length' type character.")

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
    ts_data <- unnest(ts_data, cols = 'observations')
    ts_data <- nest(ts_data, .by = 'label', .key = "observations")

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
  # add twdtw arguments to model
  model$twdtw_args <- list(time_weight = time_weight,
                           cycle_length = cycle_length,
                           time_scale = time_scale,
                           origin = NULL,
                           max_elapsed = Inf,
                           version = "f90")
  new_twdtw_args <- list(...)
  matching_twdtw_args <- intersect(names(model$twdtw_args), names(new_twdtw_args))
  model$twdtw_args[matching_twdtw_args] <- new_twdtw_args[matching_twdtw_args]

  class(model) <- "twdtw_knn1"

  return(model)

}

#' Print method for objects of class twdtw_knn1
#'
#' This method provides a structured printout of the important components
#' of a `twdtw_knn1` object.
#'
#' @param x An object of class `twdtw_knn1`.
#' @param ... ignored
#' 
#' @return Invisible `twdtw_knn1` object.
#'
#' @export
print.twdtw_knn1 <- function(x, ...) {
  cat("\nModel of class 'twdtw_knn1'\n")
  cat("-----------------------------\n")

  # Printing the call
  cat("Call:\n")
  print(x$call)

  # Printing the formula, if available
  cat("\nFormula:\n")
  print(x$formula)

  # Printing the data summary
  cat("\nData:\n")
  print(x$data)

  # Printing twdtw arguments
  cat("\nTWDTW Arguments:\n")
  pretty_arguments(x$twdtw_args)

  invisible(x) # Returns the object invisibly, so it doesn't print twice
}

#' Print Pretty Arguments
#'
#' Display a list of arguments of a given function in a human-readable format.
#'
#' @param args A list of named arguments to display.
#'
#' @return Invisible NULL. The function is mainly used for its side effect of printing.
#'
#' @examples
#' \dontrun{
#' pretty_arguments(formals(twdtw_knn1))
#' }
#'
pretty_arguments <- function(args) {

  if (is.null(args)) {
    cat("Arguments are missing.\n")
    return(invisible(NULL))
  }

  for (name in names(args)) {
    default_value <- args[[name]]

    if (is.symbol(default_value)) {
      default_value <- as.character(default_value)

    } else if (is.null(default_value)) {
      default_value <- "NULL"

    } else if (is.vector(default_value) && !is.null(names(default_value))) {
      # Handle named vectors
      values <- paste(names(default_value), default_value, sep = "=", collapse = ", ")
      default_value <- paste0("c(", values, ")")
    }
    cat(paste0(" - ", name, ": ", default_value, "\n"))
  }
}



#' Compute the Most Common Sampling Frequency across all observations
#'
#' This function calculates the most common difference between consecutive time points.
#' This can be useful for determining the aproximate sampling frequency of the time series data.
#'
#' @param x A data frame including a column called `observations`` with the time series 
#'
#' @return A difftime object representing the most common time difference between consecutive samples.
#'
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
