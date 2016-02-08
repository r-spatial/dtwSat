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
#   R Package dtwSat - 2016-02-08                             #
#                                                             #
###############################################################

#' @title Get imput time series from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves time series used in the
#' TWDTW analysis.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' 
#' @docType methods
#' @return The satellite time series, a \link[zoo]{zoo} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}.
#' 
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, 
#'              weight.fun = log_fun)
#' 
#' getTimeSeries(matches)
#' 
#' 
#' @export
setGeneric("getTimeSeries", 
           function(object){
             object@x
           }
)


