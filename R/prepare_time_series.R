#' Prepare a Time Series Tibble from a 2D stars Object with Bands and Time Attributes
#'
#' This function reshapes a data frame, which has been converted from a stars object, into a nested wide tibble format.
#' The stars object conversion often results in columns named in formats like "band.YYYY.MM.DD", "XYYYY.MM.DD.band", or "YYYY.MM.DD.band".
#'
#' @param x A data frame derived from a stars object containing time series data in wide format.
#' The column names should adhere to one of the following formats: "band.YYYY.MM.DD", "XYYYY.MM.DD.band", or "YYYY.MM.DD.band".
#'
#' @return A nested tibble in wide format. Each row of the tibble corresponds to a unique 'ts_id' that maintains the order from the original stars object.
#' The nested structure contains observations (time series) for each 'ts_id', including the 'time' of each observation, and individual bands are presented as separate columns.
#'
#'
prepare_time_series <- function(x) {

  # Remove the 'geom' column if it exists
  x$geom <- NULL
  x$x <- NULL
  x$y <- NULL
  var_names <- names(x)
  var_names <- var_names[!var_names %in% 'label']

  # Extract date and band information from the column names
  date_band <- do.call(rbind, lapply(var_names, function(name) {

    # Replace any hyphen with periods for consistent processing
    name <- gsub("-", "\\.", name)

    # Extract date and band info based on different name patterns
    if (grepl("^.+\\.\\d{4}\\.\\d{2}\\.\\d{2}$", name)) {
      date_str <- gsub("^.*?(\\d{4}\\.\\d{2}\\.\\d{2})$", "\\1", name)
      band_str <- gsub("^(.+)\\.\\d{4}\\.\\d{2}\\.\\d{2}$", "\\1", name)
    }
    else if (grepl("^X\\d{4}\\.\\d{2}\\.\\d{2}\\..+$", name)) {
      date_str <- gsub("^X(\\d{4}\\.\\d{2}\\.\\d{2})\\..+$", "\\1", name)
      band_str <- gsub("^X\\d{4}\\.\\d{2}\\.\\d{2}\\.(.+)$", "\\1", name)
    } 
    else if (grepl("^\\d{4}\\.\\d{2}\\.\\d{2}\\..+$", name)) {
      date_str <- gsub("^(\\d{4}\\.\\d{2}\\.\\d{2})\\..+$", "\\1", name)
      band_str <- gsub("^\\d{4}\\.\\d{2}\\.\\d{2}\\.(.+)$", "\\1", name)
    } 
    else {
      stop(paste("Unrecognized format in:", name))
    }

    # Convert the date string to Date format
    date <- to_date_time(gsub("\\.", "-", date_str))

    return(data.frame(time = date, band = band_str))
  }))

  # Construct tiem sereis
  ns <- nrow(x)
  x$ts_id <- 1:ns
  if (!'label' %in% names(x)) {
    x$label <- NA
  }
  x <- pivot_longer(x, !c('ts_id', 'label'), names_to = 'band_date', values_to = 'value')
  x$band <- rep(date_band$band, ns)
  x$time <- rep(date_band$time, ns)
  x$band_date <- NULL
  result_df <- pivot_wider(x, id_cols = c('ts_id', 'label', 'time'), names_from = 'band', values_from = 'value')
  result_df <- nest(result_df, .by = c('ts_id', 'label'), .key = 'observations')

  return(result_df)

}


#### TO BE REMOVED: twdtw package will export this fucntion
to_date_time <- function(x){
  if (!inherits(x, c("Date", "POSIXt"))) {
    # check if all strings in the vector include hours, minutes, and seconds
    if (all(grepl("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", x))) {
      x <- try(as.POSIXct(x), silent = TRUE)
    } else {
      x <- try(as.Date(x), silent = TRUE)
    }
    if (inherits(x, "try-error")) {
      stop("Some elements of x could not be converted to a date or datetime format")
    }
  }
  return(x)
}
