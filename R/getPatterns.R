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
#' @return a list of class \code{\link[zoo]{zoo}} 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and 
#' \code{\link[dtwSat]{getTimeSeries}}
#' 
#' @export
setGeneric("getPatterns", function(object, labels=NULL) standardGeneric("getPatterns"))

#' @rdname getPatterns
#' @aliases getPatterns-twdtwMatches
#' @examples
#' # Get patterns from objects of class twdtwTimeSeries
#' patterns = twdtwTimeSeries(timeseries=patterns.list, labels=names(patterns.list))
#' ts = twdtwTimeSeries(timeseries=example_ts.list)
#' matches = twdtwApply(x=ts, y=patterns)
#' getPatterns(matches)
#' 
#' @export
setMethod("getPatterns", c("twdtwMatches","ANY"),
          function(object, labels) getPatterns.twdtwMatches(object, labels) )

getPatterns.twdtwMatches = function(object, labels) {
  getTimeSeries(object@patterns, labels)
}

