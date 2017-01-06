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
#' @description This class stores the cross-validation. 
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwTimeSeries}}.
#' 
#' @param times Number of partitions to create.
#' 
#' @param p the percentage of data that goes to training. 
#' See \code{\link[caret]{createDataPartition}} for details.
#' 
#' @param conf.int specifies the confidence level (0-1) for interval estimation of the 
#' population mean. For more details see \code{\link[ggplot2]{mean_cl_boot}}.
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{createPatterns}} and 
#' to \code{\link[dtwSat]{twdtwApply}}.
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
#' 
#' }
NULL
setClass(
  Class = "twdtwCrossValidation",
  slots = c(partitions = "list", accuracy = "list"),
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

setGeneric("twdtwCrossValidation", 
           def = function(object, ...) standardGeneric("twdtwCrossValidation")
)

#' @inheritParams twdtwCrossValidation-class
#' @aliases twdtwCrossValidation-create
#' 
#' @describeIn twdtwCrossValidation Splits the set of time 
#' series into training and validation. The function uses stratified 
#' sampling and a simple random sampling for each stratum. For each data partition 
#' this function performs a TWDTW analysis and returns the Overall Accuracy, 
#' User's Accuracy, Produce's Accuracy, error matrix (confusion matrix), and a 
#' \code{\link[base]{data.frame}} with the classification (Predicted), the 
#' reference classes (Reference), and some TWDTW information.
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
#' cross_validation = twdtwCrossValidation(field_samples_ts, times=3, p=0.1, 
#'                           freq = 8, formula = y ~ s(x, bs="cc"), weight.fun = log_fun)
#' cross_validation
#' 
#' summary(cross_validation)
#' 
#' }
#' @export
setMethod(f = "twdtwCrossValidation", 
          definition = function(object, times, p, ...) twdtwCrossValidation.twdtwTimeSeries(object, times, p, ...))

twdtwCrossValidation.twdtwTimeSeries = function(object, times, p, ...){
  
  partitions = createDataPartition(y = labels(object), times, p, list = TRUE)
  
  res = lapply(partitions, function(I){
    training_ts = subset(object, I)
    validation_ts = subset(object, -I)
    patt = createPatterns(training_ts, ...)
    twdtw_res = twdtwApply(x = validation_ts, y = patt, n=1, ...)
    df = do.call("rbind", lapply(twdtw_res[], function(xx) xx[which.min(xx$distance),]) )
    ref = labels(twdtw_res)$timeseries
    pred = df$label
    data = data.frame(.adjustFactores(ref, pred, levels=NULL, labels=NULL), df[,!names(df)%in%"labels"])
    error.matrix = table(Predicted=data$Predicted, Reference=data$Reference)
    UA = diag(error.matrix) / rowSums(error.matrix)
    PA = diag(error.matrix) / colSums(error.matrix)
    O  = sum(diag(error.matrix)) / sum(rowSums(error.matrix))
    list(OverallAccuracy=O, UsersAccuracy=UA, ProducersAccuracy=PA, ErrorMatrix=error.matrix, data=data)
  })
  
  new("twdtwCrossValidation", partitions=partitions, accuracy=res)
  
}






