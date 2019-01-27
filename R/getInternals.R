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

setGeneric("getInternals", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getInternals"))
setGeneric("getAlignments", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getAlignments"))
setGeneric("getMatches", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getMatches"))

#' @title Get elements from twdtwMatches objects
#' @name get
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Get elements from \code{\link[dtwSat]{twdtwMatches-class}} objects. 
#' 
#' @inheritParams twdtwMatches-class
#' @param timeseries.labels a vector with labels of the time series.
#' @param patterns.labels a vector with labels of the patterns.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, and
#' \code{\link[dtwSat]{twdtwApply}}
#'
#' @examples 
#' # Getting patterns from objects of class twdtwMatches
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' mat = twdtwApply(x=ts, y=patt, weight.fun=logisticWeight(-0.1,100), keep=TRUE)
#' getPatterns(mat)
#' getTimeSeries(mat)
#' getAlignments(mat)
#' getMatches(mat)
#' getInternals(mat)
#'
#' @return a list with TWDTW results or an object \code{\link[dtwSat]{twdtwTimeSeries-class}}. 
#' 
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#'   
NULL
        
#' @aliases getAlignments
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getAlignments", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) 
            getAlignments.twdtwMatches(object, timeseries.labels, patterns.labels, attr = c("label", "from", "to", "distance", "K")) )

#' @aliases getInternals
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getInternals", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) 
            getAlignments.twdtwMatches(object, timeseries.labels, patterns.labels, attr = c("internals")) )

#' @aliases getMatches
#' @inheritParams get
#' @rdname get  
#' @export
setMethod("getMatches", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) 
            getAlignments.twdtwMatches(object, timeseries.labels, patterns.labels, attr = c("matching")) )

getAlignments.twdtwMatches = function(object, timeseries.labels, patterns.labels, attr){
  if(is.null(timeseries.labels)) timeseries.labels = labels(object@timeseries)
  if(is.null(patterns.labels)) patterns.labels = labels(object@patterns)
  res = object[timeseries.labels, patterns.labels, drop=FALSE]
  lapply(res, function(x) lapply(x, function(x) x[attr]) )
}