#'
#' Train a KNN-1 TWDTW model
#'
#' This function prepares a KNN-1 model with the Time Warp Dynamic Time Warping (TWDTW) algorithm.
#'
#' @param x A three-dimensional stars object (x, y, time) with bands as attributes.
#' @param y An sf object with the coordinates of the training points.
#' @param time_weight A numeric vector with length two (steepness and midpoint of logistic weight) or a function.
#' See details in \link[twdtw]{twdtw}.
#' @param cycle_length The length of the cycle, e.g. phenological cycles. Details in \link[twdtw]{twdtw}.
#' @param time_scale Specifies the time scale for the observations. Details in \link[twdtw]{twdtw}.
#' @param smooth_fun a function specifying how to create temporal patterns using the samples.
#' If not defined, it will keep all samples. Note that reducing the samples to patterns can significantly 
#' improve computational time of predictions. See details.
#' @param start_column Name of the column in y that indicates the start date. Default is 'start_date'.
#' @param end_column Name of the column in y that indicates the end date. Default is 'end_date'.
#' @param label_colum Name of the column in y containing land use labels. Default is 'label'.
#' @param resampling_freq The time for sampling the time series if `smooth_fun` is given.
#' If NULL, the function will infer the frequency of observations in `x`.
#' @param ... Additional arguments passed to \link[twdtw]{twdtw}.
#'
#' @details If \code{smooth_fun} not informed, the KNN-1 model will retain all training samples.
#'
#' If a custom smoothing function is passed to `smooth_fun`, the function will be used to
#' resample values of samples sharing the same label (land cover class). 
#'
#' The custom smoothing function takes two numeric vectors as arguments and returns a model:
#' \itemize{
#'   \item The first argument represents the independent variable (typically time).
#'   \item The second argument represents the dependent variable (e.g., band values) corresponding to each coordinate in the first argument.
#' }
#' See the examples section for further clarity.
#'
#' Smooting the samples can significantly reduce the processing time for prediction using `twdtw_knn1` model.
#'
#' @return A 'twdtw_knn1' model containing the trained model information and the data used.
#'
#' @examples
#' \dontrun{
#'
#' # Read training samples
#' samples_path <-
# '   system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat")
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
#' m <- twdtw_knn1(
#'    x = dc,
#'    y = samples,
#'    smooth_fun = function(x, y) gam(y ~ s(x), data = data.frame(x = x, y = y))
#'    cycle_length = 'year',
#'    time_scale = 'day',
#'    time_weight = c(steepness = 0.1, midpoint = 50))
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
#'
#' # Create a knn1-twdtw model with custom smoothing function
#'
#' m <- twdtw_knn1(x = dc,
#'  y = samples,
#'  smooth_fun = function(x, y) lm(y ~ factor(x), data = data.frame(x=x, y=y))
#'  cycle_length = 'year',
#'  time_scale = 'day',
#'  time_weight = c(steepness = 0.1, midpoint = 50))
#' 
#' plot(m)
#'
#' }
#' @export
twdtw_knn1 <- function(x, y, smooth_fun = NULL, resampling_freq = NULL,
                       time_weight, cycle_length, time_scale,
                       start_column = 'start_date', end_column = 'end_date',
                       label_colum = 'label', ...){

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

  # adjust y column names
  st_geometry(y) <- 'geom'
  y <- y[, c(start_column, end_column, label_colum, 'geom')]

  # Convert columns to date-time
  y[, start_column] <- to_date_time(y[[start_column]])
  y[, end_column] <- to_date_time(y[[end_column]])

  # Prepare time series samples from stars object
  ts_data <- st_extract(x, y)
  ts_data$label <- y[[label_colum]]
  ts_data <- prepare_time_series(as.data.frame(ts_data))
  ts_data$ts_id <- NULL

  smooth_models <- NULL
  if(!is.null(smooth_fun)) {

    # Check if smooth_fun has two or three arguments
    if(length(formals(smooth_fun)) != c(2)) {
      stop("The smooth function should have only two arguments!")
    }

    # Shift dates
    ts_data$observations <- lapply(ts_data$observations, shift_ts_dates)

    # Split data frame by label
    ts_data <- unnest(ts_data, cols = 'observations')
    ts_data <- nest(ts_data, .by = 'label', .key = "observations")

    # Apply smooth function
    smooth_models <- lapply(ts_data$observations, function(ts) {

      # Get timeline
      y_time <- ts$time
      ts$time <- NULL

      # Fit smooth model to each band
      smooth_models <- lapply(as.list(ts), function(band) {
        smooth_fun(x = as.numeric(y_time), y = band)
      })

      return(smooth_models)

    })

    names(smooth_models) <- ts_data$label

    ts_data$observations <- lapply(seq_along(ts_data$observations), function(l) {

      # Get timeline
      ts <- ts_data$observations[[l]]
      y_time <- ts$time
      ts$time <- NULL

      # Determine time for resampling time sereis
      if (is.null(resampling_freq)) {
        pred_time <- unique(y_time)
      } else {
        pred_time <- seq(min(y_time), max(y_time), by = resampling_freq)
      }

      smoothed_data <- sapply(smooth_models[[l]], function(m) {
        # Determine target class
        target_class <- class(model.frame(m)[, 2])
        # Convert pred_time based on target class
        pred_points <- if (target_class == "factor") {
          factor(as.numeric(pred_time), levels = levels(model.frame(m)[, 2]))
        } else {
          as(as.numeric(pred_time), target_class)
        }
        predict(m, newdata = data.frame(x = pred_points))
      })

      # Bind time and smoothed data into a data frame
      result_df <- data.frame(time = pred_time, smoothed_data)

      return(result_df)
    })

  }

  model <- list()
  model$call <- match.call()
  model$smooth_models <- smooth_models
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

  # Printing the smooth_fun, if available
  cat("\nRoot Mean Squared Error (RMSE) of smooth models:\n")
  if(is.null(x$smooth_models)){
    print(NULL)
  } else {
    print(sapply(x$smooth_models, function(m1) sapply(m1, function(m2){
      residuals <- try(resid(m2))
      ifelse(is.numeric(residuals), sqrt(mean(residuals^2)), NA)
    })))
  }

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
#' @noRd
#' @keywords internal
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
#' @noRd
#' @keywords internal
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
