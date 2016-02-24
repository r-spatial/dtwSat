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
#   R Package dtwSat - 2016-01-19                             #
#                                                             #
###############################################################

#' @title Join time series 
#' @name joinTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Join time series from two or more objects of class twdtwTimeSeries.
#' 
#' @param ... objects of class twdtwTimeSeries 
#'
#' @param join.labels a logical. It TRUE the function joins labels that are identical 
#' to a factor. If FALSE a different label is kept for each samples. 
#' 
#' @return an object of class \code{\link[dtwSat]{twdtwTimeSeries}}
#'
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{joinAlignments}}
#' 
#' @examples 
#' ## 
#' 
#' @export 
setGeneric(name = "joinTimeSeries",  
          def = function(...) standardGeneric("joinTimeSeries")
)

#' @rdname joinTimeSeries
#' @aliases joinTimeSeries-twdtwTimeSeries
#' @examples
#' # Join time series from objects of class twdtwTimeSeries
#' ts1 = twdtwTimeSeries(timeseries=example_ts.list[1], labels="A")
#' ts2 = twdtwTimeSeries(timeseries=example_ts.list[2], labels="B")
#' joinTimeSeries(ts1, ts2, join.labels=TRUE)
#' joinTimeSeries(ts1, ts2, join.labels=FALSE)
#' 
#' @export
setMethod("joinTimeSeries", 
          definition = function(..., join.labels=TRUE) {
            if(join.labels)
              return(joinTimeSeries.twdtwTimeSeries.join.labels(list(...)))
            joinTimeSeries.twdtwTimeSeries(list(...))
          })

joinTimeSeries.twdtwTimeSeries = function(x){
  n = length(x)
  if(n==1) return(x[[1]])
  tsn = getTimeSeries(x[[n]])
  if(length(tsn)<1) return(joinTimeSeries.twdtwTimeSeries(x[-n]))
  timeseries = c(getTimeSeries(x[[1]]), tsn)
  l1 = as.character(labels(x[[1]]))
  ln = paste(as.character(labels(x[[n]])), n, sep=".")
  labels = c(l1, ln)
  x[[1]] = twdtwTimeSeries(timeseries=timeseries, labels=labels)
  joinTimeSeries.twdtwTimeSeries(x[-n])
}

joinTimeSeries.twdtwTimeSeries.join.labels = function(x){
  n = length(x)
  if(n==1) return(x[[1]])
  tsn = getTimeSeries(x[[n]])
  if(length(tsn)<1) return(joinTimeSeries.twdtwTimeSeries.join.labels(x[-n]))
  l1 = as.character(labels(x[[1]]))
  ln = as.character(labels(x[[n]]))
  ln = ln[!ln%in%l1]
  if(length(ln)<1) return(joinTimeSeries.twdtwTimeSeries.join.labels(x[-n]))
  patterns = c(getTimeSeries(x[[1]], labels=l1), getTimeSeries(x[[n]], labels=ln))
  x[[1]] = twdtwTimeSeries(timeseries=patterns, labels=c(l1, ln))
  joinTimeSeries.twdtwTimeSeries.join.labels(x[-n])
}      


