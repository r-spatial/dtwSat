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


#' @title Get dates from year and day of the year
#' 
#' @description This function retrieves the date corresponding to the ginven 
#' year and day of the year.
#' 
#' @param year An vector with the years.
#' @param doy An vector with the day of the year.
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{Dates}} object.
#' 
#' @seealso \link[dtwSat]{getModisTimeSequence} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' year = c(2000, 2001)
#' doy = c(366, 365)
#' dates = getDatesFromDOY(year, doy)
#' dates
#'
#' @export
getDatesFromDOY = function(year, doy){
  res = as.Date(paste(as.numeric(year), as.numeric(doy)), format="%Y %j", origin="1970-01-01")
  # Correct leap years 
  I = which(diff(res)<0)
  if(length(I)>0) res[I+1] = as.Date(paste0(as.numeric(format(res[I+1], "%Y"))+1,format(res[I+1],"-%m-%d")))
  res
}

