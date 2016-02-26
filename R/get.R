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

setGeneric("getPatterns", function(object, patterns.labels=NULL) standardGeneric("getPatterns"))
setGeneric("getTimeSeries", function(object, timeseries.labels=NULL) standardGeneric("getTimeSeries"))
setGeneric("getAlignments", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getAlignments"))
setGeneric("getInternals", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getInternals"))
setGeneric("getMatches", function(object, timeseries.labels=NULL, patterns.labels=NULL) standardGeneric("getMatches"))

#' @title Get elements from twdtwMatches objects
#' @name get
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Get elements from \code{\link[dtwSat]{twdtwMatches-class}} objects. 
#' 
#' @inheritParams twdtwTimeSeries-class
#' @param timeseries.labels a vector with labels of the time series.
#' @param patterns.labels a vector with labels of the patterns.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwApply}}
#'
#' @examples 
#' # Getting patterns from objects of class twdtwMatches
#' patt = twdtwTimeSeries(patterns.list)
#' ts = twdtwTimeSeries(example_ts.list)
#' mat = twdtwApply(x=ts, y=patt, weight.fun=logisticWeight(-0.1,50), keep=TRUE)
#' getPatterns(mat)
#' getTimeSeries(mat)
#' getAlignments(mat)
#' getMatches(mat)
#' getInternals(mat)
#'
#' @return a list with TWDTW results or an object \code{\link[dtwSat]{twdtwTimeSeries-class}}. 
#'
NULL



#' @aliases getPatterns
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getPatterns", c("twdtwMatches","ANY"),
          function(object, patterns.labels) getPatterns.twdtwMatches(object, patterns.labels) )

getPatterns.twdtwMatches = function(object, patterns.labels) {
  subset(object@patterns, labels=patterns.labels)
}


#' @aliases getTimeSeries
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getTimeSeries", c("twdtwMatches","ANY"),
          function(object, timeseries.labels) getTimeSeries.twdtwMatches(object, timeseries.labels) )

getTimeSeries.twdtwMatches = function(object, timeseries.labels) {
  subset(object@timeseries, labels=timeseries.labels)
}

          
#' @aliases getAlignments
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getAlignments", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) getAlignments.twdtwMatches(object, timeseries.labels, patterns.labels) )

getAlignments.twdtwMatches = function(object, timeseries.labels, patterns.labels){
  if(is.null(timeseries.labels)) timeseries.labels = labels(object@timeseries)
  if(is.numeric(timeseries.labels)) timeseries.labels = labels(object@timeseries)[timeseries.labels]
  if(is.null(patterns.labels)) patterns.labels = labels(object@patterns)
  if(is.numeric(patterns.labels)) patterns.labels = labels(object@patterns)[patterns.labels]
  I = match(labels(object@timeseries), timeseries.labels)
  J = match(labels(object@patterns), patterns.labels)
  attr = c("label", "from", "to", "distance")
  if( all(is.na(I)) | all(is.na(J)) ) return(list())
  lapply(object@alignments[na.omit(I)], function(x) do.call("rbind", lapply(x[na.omit(J)], function(x) data.frame(x[attr]) ) ))
}


#' @aliases getInternals
#' @inheritParams get
#' @rdname get 
#' @export
setMethod("getInternals", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) getInternals.twdtwMatches(object, timeseries.labels, patterns.labels) )

getInternals.twdtwMatches = function(object, timeseries.labels, patterns.labels){
  if(is.null(timeseries.labels)) timeseries.labels = labels(object@timeseries)
  if(is.numeric(timeseries.labels)) timeseries.labels = labels(object@timeseries)[timeseries.labels]
  if(is.null(patterns.labels)) patterns.labels = labels(object@patterns)
  if(is.numeric(patterns.labels)) patterns.labels = labels(object@patterns)[patterns.labels]
  I = match(labels(object@timeseries), timeseries.labels)
  J = match(labels(object@patterns), patterns.labels)
  if( all(is.na(I)) | all(is.na(J)) ) return(list())
  lapply(object@alignments[na.omit(I)], function(x) lapply(x[na.omit(J)], function(x) x$internals) )
}


#' @aliases getMatches
#' @inheritParams get
#' @rdname get  
#' @export
setMethod("getMatches", c("twdtwMatches","ANY","ANY"),
          function(object, timeseries.labels, patterns.labels) getMatches.twdtwMatches(object, timeseries.labels, patterns.labels) )

getMatches.twdtwMatches = function(object, timeseries.labels, patterns.labels){
  if(is.null(timeseries.labels)) timeseries.labels = labels(object@timeseries)
  if(is.numeric(timeseries.labels)) timeseries.labels = labels(object@timeseries)[timeseries.labels]
  if(is.null(patterns.labels)) patterns.labels = labels(object@patterns)
  if(is.numeric(patterns.labels)) patterns.labels = labels(object@patterns)[patterns.labels]
  I = match(labels(object@timeseries), timeseries.labels)
  J = match(labels(object@patterns), patterns.labels)
  if( all(is.na(I)) | all(is.na(J)) ) return(list())
  lapply(object@alignments[na.omit(I)], function(x) lapply(x[na.omit(J)], function(x) x$matching) )
}


