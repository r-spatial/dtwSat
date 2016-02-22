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
#' @description Join time series from two or more objects of class twdtwTimeSeries 
#' 
#' @param ... objects of class twdtwTimeSeries 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @examples 
#' ## 
#' 
#' @export 
setGeneric(name = "joinTimeSeries",  
          def = function(...) standardGeneric("joinTimeSeries")
)

#' @inheritParams joinTimeSeries
#' @describeIn twdtwTimeSeries Join time series from two or more objects 
#' of class twdtwTimeSeries 
setMethod("joinTimeSeries", "twdtwTimeSeries",
          definition = function(...) joinTimeSeries.twdtwTimeSeries(list(...)))

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

          