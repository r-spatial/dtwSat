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
#   R Package dtwSat - 2016-02-22                             #
#                                                             #
###############################################################


#' @title Normalize patterns 
#' @name normalizePatterns
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Normalize the length of the temporal patterns. 
#' This function removes patterns (time series) with duplicated label. 
#' 
#' @inheritParams twdtwTimeSeries-class
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @export
setGeneric("normalizePatterns", function(object, length=NULL) standardGeneric("normalizePatterns"))

#' @inheritParams getPatterns
#' @describeIn twdtwTimeSeries Get temporal patterns from objects of class twdtwTimeSeries.
setMethod("normalizePatterns", c("twdtwTimeSeries", "numeric"),
          function(object, length=NULL) normalizePatterns.twdtwTimeSeries(object, length) )

normalizePatterns.twdtwTimeSeries = function(patterns, length){
  if(is.null(length)) length = max(nrow(patterns), na.rm=TRUE)
  lapply(patterns, function(q){
    freq = as.numeric(diff(range(index(q))))/(length-1)
    timeline = seq(min(index(q), na.rm = TRUE), max(index(q), na.rm = TRUE), by=freq)
    res = zoo(data.frame(na.spline(q, xout = timeline)), timeline)
    names(res) = names(q)
    res
  })
}

