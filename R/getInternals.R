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

#' @title Get internals from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves cost matrix, inputs, and other 
#' internal structures from \link[dtwSat]{twdtw-class} object
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} whose 
#'   elements have the internal structures used used in \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken in the alignments,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}}
#'       \cr\code{pattern}: pattern time series, 
#'       \cr\code{x}: satellite image time series,
#'       \cr\code{N}: \code{pattern} length, and 
#'       \cr\code{M}: \code{x} length.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, keep=TRUE)
#' 
#' a = getInternals(alig)
#' names(a) 
#' 
#' a = getInternals(alig, p.names="Maize")
#' names(a) 
#' 
#' a = getInternals(alig, p.names=c(1,2))
#' names(a) 
#' 
#' 
#' @export
setGeneric("getInternals", 
           function(object, ...){
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getInternals(object, p.names)
           }
)

.getInternals = function(object, p.names) {
  lapply(p.names, function(p) object@alignments[[p]]$internals)
}

