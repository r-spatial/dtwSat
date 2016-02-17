###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @title Shift dates of zoo object 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function shifts the dates of the time series to a given base year. 
#' 
#' 
#' @param x a \link[zoo]{zoo} object. 
#' @param year An integer. The base year to shift the dates of the time series to.
#' 
#' @docType methods
#' 
#' @return a \link[zoo]{zoo} object shifted to the given year.
#'
#' @examples
#' dates = seq(from = as.Date("2014-09-01"), to = as.Date("2015-03-01"), by = "15 day")
#' x = zoo(seq_along(dates), dates)
#' y = shiftDates(x, year=2005)
#' index(x)
#' index(y)
#' 
#' @export
shiftDates = function(x, year){
  dates = index(x)
  last_date = tail(dates, 1)
  shift_days = as.numeric(last_date - as.Date(paste0(year,format(last_date, "-%m-%d"))))
  d = as.numeric(dates) - shift_days
  zoo(data.frame(x), as.Date(d))
}


