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
#' @aliases twdtwTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for set of irregular time series.
#' 
#' @param ... \code{\link[dtwSat]{twdtwTimeSeries}} objects, 
#' \code{\link[zoo]{zoo}} objects or a list of \code{\link[zoo]{zoo}} objects.
#' @param labels a vector with labels of the time series. 
#' @param object an object of class twdtwTimeSeries.
#' @param x an object of class twdtwTimeSeries.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{timeseries}:}{A list of \code{\link[zoo]{zoo}} objects.}
#'  \item{\code{labels}:}{A vector of class \code{\link[base]{factor}} with time series labels.}
#' }
#'
#' @seealso   
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' \code{\link[dtwSat]{getTimeSeries}}, and
#' \code{\link[dtwSat]{twdtwApply}}
#'
#' @examples 
#' # Creating new object of class twdtwTimeSeries  
#' ptt = new("twdtwTimeSeries", timeseries = patterns.list, labels = names(patterns.list))
#' class(ptt)
#' labels(ptt)
#' levels(ptt)
#' length(ptt)
#' nrow(ptt)
#' ncol(ptt)
#' dim(ptt)
NULL
twdtwTimeSeries = setClass(
  Class = "twdtwTimeSeries",
  slots = c(timeseries = "list", labels = "factor"),
  validity = function(object){
    if(!is(object@timeseries, "list")){
      stop("[twdtwTimeSeries: validation] Invalid timeseries object, class different from list.")
    }else{}
    if(any(length(object@timeseries)>1 & !sapply(object@timeseries, is.zoo))){
      stop("[twdtwTimeSeries: validation] Invalid timeseries object, class different from list of zoo objects.")
    }else{}
    if(!is(object@labels, "factor")){
      stop("[twdtwTimeSeries: validation] Invalid labels object, class different from character.")
    }else{}
    if( length(object@labels)!=0 & length(object@labels)!=length(object@timeseries) ){
      stop("[twdtwTimeSeries: validation] Invalid labels, labels and timeseries do not have the same length.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "twdtwTimeSeries",
  definition = 
    function(.Object, timeseries, labels){
      .Object@timeseries = list()
      .Object@labels = factor(NULL)
      if(!missing(timeseries)){
        if(is(timeseries, "zoo")) timeseries = list(timeseries)
        .Object@timeseries = timeseries
        .Object@labels = factor( paste0("ts",seq_along(timeseries)) ) 
        if(!is.null(names(timeseries))) .Object@labels = factor(names(timeseries))
      }
      if(!missing(labels)){
        .Object@labels = factor(labels)
        names(.Object@timeseries) = as.character(labels)
      }
      validObject(.Object)
      return(.Object)
  }
)

setGeneric(name = "twdtwTimeSeries", 
          def = function(...) standardGeneric("twdtwTimeSeries")
)

#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Create object of class twdtwTimeSeries.
#'
#' @examples 
#' # Creating objects of class twdtwTimeSeries from zoo objects
#' ts = twdtwTimeSeries(example_ts)
#' ts 
#' 
#' # Creating objects of class twdtwTimeSeries from list of zoo objects 
#' patt = twdtwTimeSeries(patterns.list)
#' patt
#' 
#' # Joining objects of class twdtwTimeSeries 
#' tsA = twdtwTimeSeries(example_ts.list[[1]], labels = "A")
#' tsB = twdtwTimeSeries(B = example_ts.list[[2]])
#' ts = twdtwTimeSeries(tsA, tsB, C=example_ts)
#' ts
#'  
#' @export
setMethod(f = "twdtwTimeSeries", 
          definition = function(..., labels=NULL){
              timeseries = list(...)
              joint_timeseries = list()
              timeseries_class = sapply(timeseries, class)
              zoo_obj = NULL 
              list_obj = NULL 
              twdtw_obj = NULL 
              check_class = c("zoo", "list", "twdtwTimeSeries") %in% timeseries_class
              if(check_class[1]){
                  zoo_obj = timeseries[which(timeseries_class=="zoo")]
                  names(zoo_obj) = names(timeseries)[which(timeseries_class=="zoo")]
                  if(is.null(names(zoo_obj))) names(zoo_obj) = paste0("ts",seq_along(zoo_obj))
                  joint_timeseries = c(joint_timeseries, zoo_obj)
              } else {}
              if(check_class[2]){
                  list_obj = c(do.call("c", timeseries[which(timeseries_class=="list")]))
                  if(is.null(names(list_obj))) names(list_obj) = paste0("ts",seq_along(list_obj))
                  joint_timeseries = c(joint_timeseries, list_obj)
              } else {}
              if(check_class[3]){
                  twdtw_obj = do.call("c", lapply(timeseries[which(timeseries_class=="twdtwTimeSeries")], getTimeSeries))
                  names(twdtw_obj) = as.character(unlist(lapply(timeseries[which(timeseries_class=="twdtwTimeSeries")], labels)))
                  joint_timeseries = c(joint_timeseries, twdtw_obj)
              } else {}
              if(is.null(labels)) labels = names(joint_timeseries)
              new("twdtwTimeSeries", timeseries = joint_timeseries, labels = labels)
          })
          
          
          
          
          