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

#' @title Get matching points from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the matching points 
#' for each alignment between the \code{pattern} and the 
#' time series \code{x}.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves the matching points for all patterns.  
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} whose 
#'  elements have the information for each pattern. 
#'  Each pattern has then a sublist of matching points such as: 
#'       \cr\code{index1}: matching points of the pattern, and
#'       \cr\code{index2}: matching points of the time series.
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
#' getMatches(matches)
#' 
#' getMatches(matches, p.names="Maize")
#' 
#' getMatches(matches, p.names=1)
#' 
#' @export
setGeneric("getMatches", 
           function(object, p.names){
             p.names = getPatternNames(object, p.names)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getMatches(object, p.names)
           }
)

.getMatches = function(object, p.names){
  lapply(p.names, function(p) object@alignments[[p]]$matching)
}

