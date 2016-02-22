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


#' @title class "twdtwTimeSeries"
#' @name twdtwTimeSeries-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for set of irregular time series.
#' 
#' @param ... objects of class \code{\link[zoo]{zoo}}.
#' @param timeseries a list of \code{\link[zoo]{zoo}} objects.
#' @param labels a vector with labels of the time series. 
#' @param object an object of class twdtwTimeSeries.
#' @param x an object of class twdtwTimeSeries.
#' @param length a number. If not informed the function will use the 
#' length of the longest pattern. 
#'
#' @section Slots :
#' \describe{
#'  \item{\code{timeseries}:}{A list of \code{\link[zoo]{zoo}} objects.}
#'  \item{\code{labels}:}{A vector of class \code{\link[base]{factor}} with time series labels.}
#' }
#' 
#' @examples 
#' # Creating new object of class twdtwTimeSeries  
#' new("twdtwTimeSeries")
#' 
#' sapply(example_ts.list, class)
#' new("twdtwTimeSeries", timeseries = example_ts.list)
#' new("twdtwTimeSeries", timeseries = example_ts.list, labels = c(1,2))
#' 
NULL
twdtwTimeSeries = setClass(
  Class = "twdtwTimeSeries",
  slots = c(timeseries = "list", labels = "factor"),
  validity = function(object){
    if(!is(object@timeseries, "list")){
      stop("[twdtwTimeSeries: validation] Invalid timeseries object, class different from list.")
    }else{}
    if(any(!sapply(object@timeseries, is.zoo))){
      stop("[twdtwTimeSeries: validation] Invalid timeseries object, class different from list of zoo objects.")
    }else{}
    if(!is(object@labels, "factor")){
      stop("[twdtwTimeSeries: validation] Invalid labels object, class different from character.")
    }else{}
    if( length(object@labels)!=0 & length(object@labels)!=length(object@timeseries) ){
      stop("[twdtwTimeSeries: validation] Invalid length, labels and timeseries do not have the same length.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "twdtwTimeSeries",
  definition = 
    function(.Object, timeseries, labels){
      .Object@timeseries = list(zoo(NULL))
      .Object@labels = factor(NULL)
      if(!missing(timeseries)){
        if(is(timeseries, "zoo")) timeseries = list(timeseries)
        .Object@timeseries = timeseries
        .Object@labels = factor(paste0("ts",seq_along(.Object@timeseries)))
      }
      if(!missing(labels))
        .Object@labels = factor(labels)
      validObject(.Object)
      return(.Object)
  }
)

#' @title Create twdtwTimeSeries object 
#' @name twdtwTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Create object of class twdtwTimeSeries.
#' 
#' @inheritParams twdtwTimeSeries-class
#' 
#' @examples 
#' # Creating objects of class twdtwTimeSeries 
#' class(example_ts)
#' twdtwTimeSeries(example_ts)
#' sapply(example_ts.list, class)
#' twdtwTimeSeries(timeseries = example_ts.list)
#' twdtwTimeSeries(timeseries = example_ts.list, labels = c(1,2))
#' twdtwTimeSeries(example_ts.list[[1]], example_ts.list[[2]])
#' 
#' @export
setGeneric(name = "twdtwTimeSeries",  
          def = function(..., timeseries, labels) standardGeneric("twdtwTimeSeries")
)

#' @inheritParams twdtwTimeSeries
#' @describeIn twdtwTimeSeries Create object of class twdtwTimeSeries.
setMethod(f = "twdtwTimeSeries",  
          definition = function(..., timeseries = list(...), labels)
             new("twdtwTimeSeries", timeseries = timeseries, labels = labels)
          )

