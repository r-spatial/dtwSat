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
#' for each alignment between the \code{pattern} and 
#' \code{x}
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the temporal pattern and the time series. 
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the pattern, and
#'       \cr\code{index2}: matching points of the time series.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun)
#' 
#' getMatches(alig)
#' 
#' getMatches(alig, p.names="Maize")
#' 
#' getMatches(alig, p.names=1)
#' 
#' @export
setGeneric("getMatches", 
           function(object, ...){
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getMatches(object, p.names)
           }
)

.getMatches = function(object, p.names){
  lapply(p.names, function(p) object@alignments[[p]]$matching)
}

