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
#   R Package dtwSat - 2017-01-18                             #
#                                                             #
###############################################################


#' @title class "twdtwAssessment" 
#' @name twdtwAssessment-class
#' @aliases twdtwAssessment
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This class stores the map assessment metrics.  
#' 
#' @param object an object of class twdtwAssessment.
#' 
#' @seealso \code{\link[dtwSat]{twdtwClassify}},  
#' \code{\link[dtwSat]{twdtwAssess}}, and
#' \code{\link[dtwSat]{twdtwXtable}}.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{accuracySummary}:}{Overall Accuracy, User's Accuracy, Produce's Accuracy, 
#'  Error Matrix (confusion matrix), and Estimated Area, considering all time periods.}
#'  \item{\code{accuracyByPeriod}:}{Overall Accuracy, User's Accuracy, Produce's Accuracy, 
#'  Error Matrix (confusion matrix), and Estimated Area, for each time periods independently 
#'  from each other.}
#'  \item{\code{data}:}{A \code{\link[base]{data.frame}} with period (from - to), reference labels, 
#'  predicted labels, and other TWDTW information.}
#' }
#' 
#' @details
#' If the twdtwRaster is unprojected (longitude/latitude) the estimated area is sum of the approximate 
#' surface area in km2 of each cell (pixel). If the twdtwRaster is projected the estimated area is calculated 
#' using the the pixel resolution in the map unit.
#'
NULL
setClass(
  Class = "twdtwAssessment",
  slots = c(accuracySummary = "list", accuracyByPeriod = "list", data = "data.frame"),
  validity = function(object){
    if(!is(object@accuracySummary, "list")){
      stop("[twdtwAssessment: validation] Invalid partitions, class different from list.")
    }else{}
    if(!is(object@accuracyByPeriod, "list")){
      stop("[twdtwAssessment: validation] Invalid accuracy, class different from list.")
    }else{}
    if(!is(object@data, "data.frame")){
      stop("[twdtwAssessment: validation] Invalid accuracy, class different from data.frame.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
          signature = "twdtwAssessment",
          definition = 
            function(.Object, accuracySummary, accuracyByPeriod, data){
              .Object@accuracySummary = list(OverallAccuracy=NULL, UsersAccuracy=NULL, ProducersAccuracy=NULL, ErrorMatrix=table(NULL))
              .Object@accuracyByPeriod = list(list(OverallAccuracy=NULL, UsersAccuracy=NULL, ProducersAccuracy=NULL, 
                                                   ErrorMatrix=table(NULL)))
              .Object@data = data.frame(Period=NULL, from=NULL, to=NULL, Distance=NULL, Predicted=NULL, Reference=NULL)
              if(!missing(accuracySummary))
                .Object@accuracySummary = accuracySummary
              if(!missing(accuracyByPeriod))
                .Object@accuracyByPeriod = accuracyByPeriod
              if(!missing(data))
                .Object@data = data
              validObject(.Object)
              return(.Object)
            }
)


