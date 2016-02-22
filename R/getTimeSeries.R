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

#' @title Get time series 
#' @name getTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Generic method to get time series from objects of class twdtw*.
#' 
#' @param object an twdtw* object.
#' @param labels a vector with the time series labels. If not informed the 
#' function retrieves all time series. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @export
setGeneric("getTimeSeries", function(object, labels=NULL) standardGeneric("getTimeSeries"))

#' @inheritParams getTimeSeries
#' @describeIn twdtwTimeSeries Get time series from objects of class twdtwTimeSeries 
setMethod("getTimeSeries", c("twdtwTimeSeries", "ANY"),
          function(object, labels) getTimeSeries.twdtwTimeSeries(object, labels) )

#' @inheritParams getTimeSeries
#' @describeIn twdtwMatches Get time series from objects of class twdtwMatches 
setMethod("getTimeSeries", c("twdtwMatches", "ANY"),
          function(object, labels) getTimeSeries.twdtwMatches(object, labels) )

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
