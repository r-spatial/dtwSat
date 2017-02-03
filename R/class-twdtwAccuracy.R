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
#' @description This class stores the map assessment.  
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwRaster}} resulting from 
#' the classification, i.e. \code{\link[dtwSat]{twdtwClassify}}.
#' 
#' @param y a \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval.
#' If missing ''from'' and/or ''to'', they are set to the time range of the 
#' \code{object}. 
#' 
#' @param id.labels a numeric or character with an column name from \code{y} to 
#' be used as samples labels. Optional.
#' 
#' @param labels character vector with time series labels. For signature 
#' \code{\link[dtwSat]{twdtwRaster}} this argument can be used to set the 
#' labels for each sample in \code{y}, or it can be combined with \code{id.labels} 
#' to select samples with a specific label.
#' 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}.
#' 
#' @param conf.int specifies the confidence level (0-1).
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}}, and 
#' \code{\link[dtwSat]{twdtwClassify}}.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{accuracySummary}:}{Overall Accuracy, User's Accuracy, Produce's Accuracy, 
#'  and Error Matrix (confusion matrix) considering all time periods.}
#'  \item{\code{accuracyByPeriod}:}{Overall Accuracy, User's Accuracy, Produce's Accuracy, 
#'  and Error Matrix (confusion matrix) for each time periods independently from each other.}
#'  \item{\code{data}:}{A \code{\link[base]{data.frame}} with period (from - to), reference labels, 
#'  predicted labels, and other TWDTW information.}
#' }
#'
#' @examples
#' \dontrun{
#' 
#' }
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


