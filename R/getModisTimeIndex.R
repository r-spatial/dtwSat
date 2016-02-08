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


#' @title Get Modis time index from date
#' 
#' @description This function retrieves the nearest time 
#' index to a date
#' 
#' @param date A \code{\link[base]{Dates}} object
#' @param frequency An integer with the frequency in days. Default is 16 days
#' 
#' @docType methods 
#' 
#' @return An \code{\link[base]{integer}} object.
#' 
#' @seealso 
#' \link[dtwSat]{getModisTimeSequence}, and
#' \link[dtwSat]{getDatesFromDOY}.
#' 
#' @examples
#' i = getModisTimeIndex(date=as.Date("2000-01-01"))
#' i
#'
#' @export
getModisTimeIndex = function(date, frequency=16){
  dates = getModisTimeSequence(frequency=frequency)
  res = which.min(abs(dates-date))
  res
}

