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


#' @title Create time sequence
#' 
#' @description This function creates a sequence of dates for 
#' each year. The sequences start on January 1st of each year.
#' 
#' @param year A vector with the years. Default 
#' is form 2000 to the system time year \code{format(Sys.time(), ''\%Y'')}
#' @param frequency An integer with the frequency in days. Default is 16 days
#' @docType methods
#' 
#' @return A vector of \code{\link[base]{Dates}} 
#' 
#' @seealso \link[dtwSat]{getDatesFromDOY} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' dates = getModisTimeSequence()
#' dates
#' 
#' @export
getModisTimeSequence = function(year=2000:format(Sys.time(), "%Y"), frequency=16){
  res = unlist(lapply(year, function(y){
    days = seq(from = as.Date(paste0(y,"-01-01")), to = as.Date(paste0(y,"-12-31")), by = frequency)
  }))
  res = as.Date(res, origin="1970-01-01")
  res
}


