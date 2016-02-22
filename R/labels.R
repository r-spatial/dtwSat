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


#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Labels of the time series in object of class twdtwTimeSeries.
#' @export
setMethod("labels", signature = signature(object="twdtwTimeSeries"),
          definition = function(object) object@labels)

#' @inheritParams twdtwMatches-class
#' @describeIn twdtwMatches Labels of the time series in object of class twdtwMatches.
#' @export
setMethod("labels", 
          signature = signature(object="twdtwMatches"),
          definition = function(object){
            list(timeseries = labels(object@timeseries), 
                 patterns = labels(object@patterns))
          }
)

