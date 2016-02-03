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
#   R Package dtwSat - 2016-30-01                             #
#                                                             #
###############################################################

#' @title Get the number of matches from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the number of matches  
#' from a \link[dtwSat]{twdtw-class} object.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves the total number of matches considering all patterns.  
#' 
#' @docType methods
#' @return A \link[base]{numeric} 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun)
#' 
#' nmatches(matches)
#' 
#' nmatches(matches, p.names="Soybean")
#' 
#' nmatches(matches, p.names=c(2,3))
#' 
#' @export
setGeneric("nmatches", 
           function(object, p.names) {
               x = getAlignments(object, p.names)
               nrow(x)
           }
)


