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


setGeneric("is.twdtwTimeSeries", 
           function(x, ...) standardGeneric("is.twdtwTimeSeries"))

setGeneric("is.twdtwMatches", 
           function(x, ...) standardGeneric("is.twdtwMatches"))

#' @aliases is.twdtwTimeSeries
#' @inheritParams twdtwTimeSeries-class
#' 
#' @examples 
#' # Checking if object belongs to the class twdtwTimeSeries
#' ex_ts = twdtwTimeSeries(timeseries = example_ts.list)
#' is.twdtwTimeSeries(ex_ts)
#' 
#' @describeIn twdtwTimeSeries Check if the object belongs to the class twdtwTimeSeries
#' @export
setMethod("is.twdtwTimeSeries", "twdtwTimeSeries", 
          function(x) is(x, "twdtwTimeSeries"))

#' @aliases is.twdtwMatches
#' @inheritParams twdtwMatches-class
#' @examples 
#' # Checking if object belongs to the class twdtwMatches
#' matches = new("twdtwMatches")
#' is.twdtwMatches(matches)
#' 
#' @describeIn twdtwMatches Check if the object belongs to the class twdtwMatches
#' @export
setMethod("is.twdtwMatches", "twdtwMatches", 
          function(x) is(x, "twdtwMatches"))

