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
#   R Package dtwSat - 2016-02-18                             #
#                                                             #
###############################################################

#' @include getTimeSeries.R
#' @title Get patterns 
#' @name getPatterns
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Generic method to get temporal patterns from objects of class twdtwMatches.
#' 
#' @param object an twdtwMatches object.
#' @param labels a vector with the patterns labels. If not informed the 
#' function retrieves all patterns. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @export
setGeneric("getPatterns", function(object, labels) standardGeneric("getPatterns"))

#' @inheritParams getPatterns
#' @describeIn twdtwMatches Get temporal patterns from objects of class twdtwMatches.
setMethod("getPatterns", c("twdtwTimeSeries", "ANY"),
          function(object, labels=NULL) getPatterns.twdtwMatches(object, labels) )

getPatterns.twdtwMatches = function(object, labels) {
  getTimeSeries(object@patterns, labels)
}

