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

#' @name getTimeSeries-method 
#' @aliases getTimeSeries
#' 
#' @title Get time series 
#' 
#' @description Generic method to get time series from objects of class twdtw*.
#' 
#' @inheritParams twdtwTimeSeries-class
#' @param object an twdtw* object.
#' @param ... other arguments to be passed to methods.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
setGeneric("getTimeSeries", function(object, ...) standardGeneric("getTimeSeries"))


#' @name getPatterns-method 
#' @aliases getPatterns
#' @include methods-getTimeSeries.R
#' 
#' @title Get time series 
#' 
#' @description Generic method to get temporal patterns from objects of class twdtw*.
#' 
#' @inheritParams twdtwTimeSeries-class
#' @param object an twdtw* object.
#' @param ... other arguments to be passed to methods.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#'
setGeneric("getPatterns", function(object, ...) standardGeneric("getPatterns"))


# Get patterns from object of class twdtwMatches by labels 
getPatterns.twdtwMatches = function(object, labels) {
  getTimeSeries(object@patterns, labels)
}

# Get time series from object of class twdtwTimeSeries by labels 
getTimeSeries.twdtwTimeSeries = function(object, labels){
  if(is.null(labels)) labels = labels(object)
  if(is.numeric(labels)) labels = labels(object)[labels]
  I = match(object@labels, labels)
  object@timeseries[na.omit(I)]
}

# Get time series from object of class twdtwMatches by labels 
getTimeSeries.twdtwMatches = function(object, labels) {
  getTimeSeries.twdtwTimeSeries(object@timeseries, labels)
}


#' @inheritParams getTimeSeries-method
#' 
#' @examples 
#' # Getting time series from objects of class twdtwTimeSeries
#' ex_ts = twdtwTimeSeries(timeseries=example_ts.list)
#' getTimeSeries(ex_ts, 1)
#' 
#' @describeIn twdtwTimeSeries Get time series from objects of class twdtwTimeSeries
#' @export
setMethod("getTimeSeries", "twdtwTimeSeries", 
          function(object, labels=NULL) getTimeSeries.twdtwTimeSeries(object, labels) )


#' @inheritParams twdtwMatches-class
#' 
#' @examples 
#' ## Getting time series from objects of class twdtwMatches 
#' matches = new("twdtwMatches")
#' getTimeSeries(matches)
#' 
#' @describeIn twdtwMatches Get time series from objects of class twdtwMatches
#' @export
setMethod("getTimeSeries", "twdtwMatches", 
          function(object, labels=NULL) getTimeSeries.twdtwMatches(object, labels) )


#' @inheritParams twdtwMatches-class
#' 
#' @examples 
#' ## Getting patterns from objects of class twdtwMatches 
#' matches = new("twdtwMatches")
#' getPatterns(matches)
#' 
#' @describeIn twdtwMatches Get patterns from objects of class twdtwMatches
#' @export
setMethod("getPatterns", "twdtwMatches", 
          function(object, labels=NULL) getPatterns.twdtwMatches(object, labels) )

