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
#' @examples 
#' # Getting labels of time series in objects of class twdtwTimeSeries
#' ex_ts = twdtwTimeSeries(timeseries = example_ts.list)
#' length(ex_ts)
#' 
#' @describeIn twdtwTimeSeries Labels of the time series in object of class twdtwTimeSeries
#' @export
setMethod("labels", 
          signature = signature(object="twdtwTimeSeries"),
          definition = function(object){
            object@labels
          }
)

#' @inheritParams twdtwMatches-class
#' @examples 
#' # Getting labels of time series in objects of class twdtwMatches
#' matches = new("twdtwMatches")
#' labels(matches)
#' 
#' @describeIn twdtwMatches Labels of the time series in object of class twdtwMatches
#' @export
setMethod("labels", 
          signature = signature(object="twdtwMatches"),
          definition = function(object){
            list(timeseries = labels(object@timeseries), 
                 patterns = labels(object@patterns))
          }
)

