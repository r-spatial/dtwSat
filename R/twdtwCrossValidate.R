
setGeneric("twdtwCrossValidate", 
           def = function(object, ...) standardGeneric("twdtwCrossValidate")
)

#' @title Cross Validate temporal patterns  
#' @name twdtwCrossValidate
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Splits the set of time series into training and validation and 
#' compute accuracy metrics. The function uses stratified sampling and a simple 
#' random sampling for each stratum. For each data partition this function 
#' performs a TWDTW analysis and returns the Overall Accuracy, User's Accuracy, 
#' Produce's Accuracy, error matrix (confusion matrix), and a \code{\link[base]{data.frame}} 
#' with the classification (Predicted), the reference classes (Reference), 
#' and the results of the TWDTW analysis.
#'
#' @param object an object of class \code{\link[dtwSat]{twdtwTimeSeries}}.
#' 
#' @param times Number of partitions to create.
#' 
#' @param p the percentage of data that goes to training. 
#' See \code{\link[caret]{createDataPartition}} for details.
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{createPatterns}} and 
#' to \code{\link[dtwSat]{twdtwApply}}.
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
#' twdtwXtable(cross_validation)
#' 
#' twdtwXtable(cross_validation, show.overall=FALSE)
#' 
#' }
NULL

#' @aliases twdtwCrossValidate-twdtwTimeSeries
#' @inheritParams twdtwCrossValidate
#' @rdname twdtwCrossValidate 
#' @export
setMethod(f = "twdtwCrossValidate", signature = "twdtwTimeSeries",
          definition = function(object, times, p, ...) twdtwCrossValidate.twdtwTimeSeries(object, times, p, ...))

twdtwCrossValidate.twdtwTimeSeries = function(object, times, p, ...){
  
  partitions = createDataPartition(y = labels(object), times, p, list = TRUE)
  
  res = lapply(partitions, function(I){
    training_ts = subset(object, I)
    validation_ts = subset(object, -I)
    patt = createPatterns(training_ts, ...)
    twdtw_res = twdtwApply(x = validation_ts, y = patt, n=1, ...)
    df = do.call("rbind", lapply(twdtw_res[], function(xx) {
      i = which.min(xx$distance)
      if(length(i)<1)
        return(data.frame(Alig.N=NA, from=NA, to=NA, distance=NA, label = "Unclassified"))
      xx[i,]
    }))
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






