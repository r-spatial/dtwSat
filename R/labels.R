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
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("labels", signature = signature(object="twdtwTimeSeries"),
          definition = function(object) object@labels)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod("labels", signature = signature(object="twdtwRaster"),
          definition = function(object) labels(object@labels))
          
#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod("labels", 
          signature = signature(object="twdtwMatches"),
          definition = function(object){
            list(timeseries = labels(object@timeseries), 
                 patterns = labels(object@patterns))
          }
)

