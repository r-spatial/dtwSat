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
#   R Package dtwSat - 2016-11-27                             #
#                                                             #
###############################################################


#' @title class "twdtwCrossValidation" 
#' @name twdtwCrossValidation-class
#' @aliases twdtwCrossValidation
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This class stores the results of the cross-validation. 
#' 
#' @param object an object of class twdtwCrossValidation. 
#' 
#' @param conf.int specifies the confidence level (0-1) for interval estimation of the 
#' population mean. For more details see \code{\link[ggplot2]{mean_cl_boot}}.
#' 
#' @param ... Other arguments. Not used. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}},
#' \code{\link[dtwSat]{createPatterns}}, and 
#' \code{\link[dtwSat]{twdtwApply}}.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{partitions}:}{A list with the indices of time series used for training.}
#'  \item{\code{accuracy}:}{A list with the accuracy and other TWDTW information for each 
#'  data partitions.}
#' }
#'
#' @examples 
#' \dontrun{
#' # Data folder 
#' data_folder = system.file("lucc_MT/data", package = "dtwSat")
#' 
#' # Read dates 
#' dates = scan(paste(data_folder,"timeline", sep = "/"), what = "dates")
#' 
#' # Read raster time series 
#' evi = brick(paste(data_folder,"evi.tif", sep = "/"))
#' raster_timeseries = twdtwRaster(evi, timeline = dates)
#' 
#' # Read field samples 
#' field_samples = read.csv(paste(data_folder,"samples.csv", sep = "/")) 
#' table(field_samples[["label"]])
#' 
#' # Read field samples projection 
#' proj_str = scan(paste(data_folder,"samples_projection", sep = "/"), 
#'      what = "character")
#' 
#' # Get sample time series from raster time series 
#' field_samples_ts = getTimeSeries(raster_timeseries, 
#'      y = field_samples, proj4string = proj_str)
#' field_samples_ts
#' 
#' # Run cross validation
#' set.seed(1)
#' # Define TWDTW weight function 
#' log_fun = logisticWeight(alpha=-0.1, beta=50) 
#' cross_validation = twdtwCrossValidate(field_samples_ts, times=3, p=0.1, 
#'                           freq = 8, formula = y ~ s(x, bs="cc"), weight.fun = log_fun)
#' cross_validation
#' 
#' summary(cross_validation)
#' 
#' plot(cross_validation)
#' 
#' }
NULL
setClass(
  Class = "twdtwCrossValidation",
  slots = c(partitions = "list", accuracy = "list"),
  validity = function(object){
    if(!is(object@partitions, "list")){
      stop("[twdtwCrossValidation: validation] Invalid partitions, class different from list.")
    }else{}
    if(!is(object@accuracy, "list")){
      stop("[twdtwCrossValidation: validation] Invalid accuracy, class different from list.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
          signature = "twdtwCrossValidation",
          definition = 
            function(.Object, partitions, accuracy){
              .Object@partitions = list(Resample1=NULL)
              .Object@accuracy = list(OverallAccuracy=NULL, UsersAccuracy=NULL, ProducersAccuracy=NULL, 
                                      ErrorMatrix=table(NULL), data=data.frame(NULL))
              if(!missing(partitions))
                .Object@partitions = partitions
              if(!missing(accuracy))
                .Object@accuracy = accuracy
              validObject(.Object)
              return(.Object)
            }
)

