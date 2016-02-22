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

#' @title Join patterns 
#' @name joinPatterns
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Join patterns (time series) from two or more objects of class twdtwTimeSeries. 
#' This function removes patterns (time series) with duplicated label. 
#' 
#' @param ... objects of class twdtwTimeSeries 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}
#' 
#' @examples 
#' # Joining temporal patterns 
#' 
#' 
#' @export
setGeneric("joinPatterns", function(...) {
  standardGeneric("joinPatterns")
})

#' @inheritParams joinPatterns
#' @describeIn twdtwTimeSeries Join patterns (time series) from two or more 
#' objects of class twdtwTimeSeries. This function removes patterns (time series) 
#' with duplicated label. 
setMethod("joinPatterns", "twdtwTimeSeries",
  function(...) joinPatterns.twdtwTimeSeries(list(...))
)

joinPatterns.twdtwTimeSeries = function(x){
  n = length(x)
  if(n==1) return(x[[1]])
  tsn = getTimeSeries(x[[n]])
  if(length(tsn)<1) return(joinPatterns.twdtwTimeSeries(x[-n]))
  l1 = as.character(labels(x[[1]]))
  ln = as.character(labels(x[[n]]))
  ln = ln[!ln%in%l1]
  if(length(ln)<1) return(joinPatterns.twdtwTimeSeries(x[-n]))
  patterns = c(getTimeSeries(x[[1]], labels=l1), getTimeSeries(x[[n]], labels=ln))
  x[[1]] = twdtwTimeSeries(timeseries=patterns, labels=c(l1, ln))
  joinPatterns.twdtwTimeSeries(x[-n])
}