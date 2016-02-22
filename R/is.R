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
           function(x) standardGeneric("is.twdtwTimeSeries"))

setGeneric("is.twdtwMatches", 
           function(x) standardGeneric("is.twdtwMatches"))

#' @aliases is.twdtwTimeSeries
#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Check if the object belongs to the class twdtwTimeSeries.
#' @export
setMethod("is.twdtwTimeSeries", "ANY", 
          function(x) is(x, "twdtwTimeSeries"))

#' @aliases is.twdtwMatches
#' @inheritParams twdtwMatches-class
#' @describeIn twdtwMatches Check if the object belongs to the class twdtwMatches.
#' @export
setMethod("is.twdtwMatches", "ANY", 
          function(x) is(x, "twdtwMatches"))

