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
#   R Package dtwSat - 2016-16-01                             #
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
#' @return a \link[zoo]{zoo} object shifted to the given year 
#'
#' @examples
#' dates = seq(from = as.Date("2014-09-01"), to = as.Date("2015-03-01"), by = "15 day")
#' x = zoo(seq_along(dates), dates)
#' y = shiftDates(x, year=2005)
#' index(x)
#' index(y)
#' 
#' @export
#' 
shiftDates = function(x, year){
  doy = as.numeric(format(index(x), "%j"))
  d = .shiftDates(year, doy)
  index(x) = d
  x
}

.shiftDates = function(year, doy){
  dates = getDatesFromDOY(year, doy)
  diffdate = diff(dates)
  I = which(diffdate<0)
  if(length(I)<1) return(dates)
  dates = c(.shiftDates(year-1, doy[1:(I[1])]), dates[(I[1]+1):length(doy)])
  dates
}

