#' Shift Dates to Start on a Specified Origin Year
#'
#' Shifts a vector of dates to start on the same day-of-year in a specified origin year
#' while preserving the relative difference in days among the observations. 
#' This way the temporal pattern (e.g., seasonality) inherent to the original dates 
#' will also be preserved in the shifted dates.
#'
#' The primary goal of this function is to align a sequence of dates based on the day-of-year
#' in a desired origin year. This can be particularly useful for comparing or visualizing
#' two or more time series with different absolute dates but aiming to align them based on
#' the day-of-year or another relative metric.
#'
#' @param x A vector of date strings or Date objects representing the sequence to shift.
#' @param origin A date string or Date object specifying the desired origin year for the shifted dates.
#' Default is "1970-01-01".
#'
#' @return A vector of Date objects with the shifted dates starting on the same day-of-year in the specified origin year.
#'
#' @examples
#' 
#' x <- c("2011-09-14", "2011-09-30", "2011-10-16", "2011-11-01")
#' 
#' shift_dates(x)
#' 
#' @export 
shift_dates <- function(x, origin = "1970-01-01") {
  
  # Convert the input dates to Date objects
  x <- as.Date(x)
  
  # Extract the day-of-year from the first date in x
  doy <- as.numeric(format(x[1], "%j"))
  
  # Compute the new starting date in the origin year
  new_start <- as.Date(paste0(format(as.Date(origin), "%Y"), sprintf("-%03d", doy)), format="%Y-%j")
  
  # Determine the difference in days
  day_diff <- as.numeric(x[1] - new_start)
  
  # Adjust each date in x by this difference
  shifted_dates <- x - day_diff
  
  return(shifted_dates)
}


shift_ts_dates <- function(x) {
  x$time <- shift_dates(x$time)
  return(x)
}
