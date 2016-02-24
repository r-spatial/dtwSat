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


#' @title Resample time series 
#' @name resampleTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description resample time series in the same object to have the same 
#' the length. 
#' 
#' @inheritParams twdtwTimeSeries-class
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwApply}}
#'
#' @return An object of class \code{\link[dtwSat]{twdtwTimeSeries}} whose 
#' time series have the same number of samples (points). 
#'
#' @export
setGeneric("resampleTimeSeries", function(object, length=NULL) standardGeneric("resampleTimeSeries"))

#' @rdname resampleTimeSeries
#' @aliases resampleTimeSeries-twdtwMatches
#' @examples
#' # Getting time series from objects of class twdtwTimeSeries
#' patterns = twdtwTimeSeries(timeseries=patterns.list, labels=names(patterns.list))
#' npatterns = resampleTimeSeries(patterns, length=46)
#' nrow(patterns)
#' nrow(npatterns)
#' 
#' @export
setMethod("resampleTimeSeries", "twdtwTimeSeries",
          function(object, length) {
              if(is.null(length)) length = max(nrow(object), na.rm=TRUE)
              do.call("joinTimeSeries", lapply(as.list(object), resampleTimeSeries.twdtwTimeSeries, length=length)) 
          })
            
resampleTimeSeries.twdtwTimeSeries = function(x, length){
  labels = as.character(labels(x))
  x = x[[1]]
  dates = index(x)
  freq = as.numeric(diff(range(dates)))/(length-1)
  timeline = seq(min(dates, na.rm = TRUE), max(dates, na.rm = TRUE), by=freq)
  res = zoo(data.frame(na.spline(x, xout = timeline)), timeline)
  new("twdtwTimeSeries", timeseries=res, labels=labels)
}






