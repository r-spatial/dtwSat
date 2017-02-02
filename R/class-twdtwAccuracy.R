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
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}.
#' 
#' @param conf.int specifies the confidence level.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}}, and 
#' \code{\link[dtwSat]{twdtwClassify}}.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{accuracy}:}{A list with the accuracy for each classified time period.}
#'  \item{\code{data}:}{A \code{\link[base]{data.frame}} with reference labels, predicted labels, 
#'  and other TWDTW information.}
#' }
#'
#' @examples
#' \dontrun{
#' 
#' }
NULL
setClass(
  Class = "twdtwAssessment",
  slots = c(accuracy = "list", data = "list"),
  validity = function(object){
    if(!is(object@partitions, "list")){
      stop("[twdtwTimeSeries: validation] Invalid partitions, class different from list.")
    }else{}
    if(!is(object@accuracy, "list")){
      stop("[twdtwTimeSeries: validation] Invalid accuracy, class different from list.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
          signature = "twdtwAssessment",
          definition = 
            function(.Object, partitions, accuracy){
              .Object@partitions = list(Resample1=NULL)
              .Object@accuracy = list(OverallAccuracy=NULL, UsersAccuracy=NULL, ProducersAccuracy=NULL, 
                                      error.matrix=table(NULL), data=data.frame(NULL))
              if(!missing(partitions))
                .Object@partitions = partitions
              if(!missing(accuracy))
                .Object@accuracy = accuracy
              validObject(.Object)
              return(.Object)
            }
)


