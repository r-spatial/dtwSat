#' @title Get dates from year and day of the year
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the date corresponding to the given 
#' year and day of the year.
#' 
#' @param year A vector with the years.
#' @param doy A vector with the day of the year. 
#' It must have the same length as \code{year}.
#' 
#' @return A \code{\link[base]{Dates}} object.
#' 
#'   
#' @export
get_dates_from_doy = function(year, doy){
  res = as.Date(paste(as.numeric(year), as.numeric(doy)), format="%Y %j", origin="1970-01-01")
  I = which(diff(res)<0)+1
  res[I] = as.Date(paste0(as.numeric(format(res[I],"%Y"))+1, format(res[I], "-%m-%d")))
  res
}



#' @title Shift dates 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function shifts the dates of the time series to a 
#' given base year. 
#' 
#' @param object \code{\link[base]{data.frame}} objects.
#' 
#' @param year the base year to shift the time series to. 
#' 
#' @return An object of the same class as the input \code{object}. 
#'
#' @export
shift_dates = function(x, year){
  labels = as.character(labels(x))
  x = x[[1]]
  dates = index(x)
  last_date = tail(dates, 1)
  shift_days = as.numeric(last_date - as.Date(paste0(year,format(last_date, "-%m-%d"))))
  d = as.numeric(dates) - shift_days
  new("twdtwTimeSeries", timeseries=zoo(data.frame(x), as.Date(d)), labels=labels)
}
