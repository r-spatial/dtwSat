#' @title Cross-validation
#' @name Cross-validation
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This functions create data partitions and compute Cross-validation metrics. 
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwTimeSeries}} or 
#' \code{\link[dtwSat]{twdtwMatches}}.
#' 
#' @param times Number of partitions to create.
#' 
#' @param p the percentage of data that goes to training. 
#' See \code{\link[caret]{createDataPartition}} for details.
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{createPatterns}}.
#' 
#' @param matrix logical. If TRUE retrieves the confusion matrix. 
#' FALSE retrieves User's Accuracy (UA) and Producer's Accuracy (PA). 
#' Dafault is FALSE. 
#' 
#' @details 
#' \describe{
#'  \item{\code{splitDataset}:}{This function splits the a set of time 
#'        series into training and validation. The function uses stratified 
#'        sampling and a simple random sampling for each stratum. Each data partition 
#'        returned by this function has the temporal patterns and a set of time series for 
#'        validation.}
#'  \item{\code{twdtwCrossValidation}:}{The function \code{splitDataset} performs the Cross-validation of 
#'        the classification based on the labels of the classified time series 
#'        (Reference) and the labels of the classification (Predicted). This function
#'        returns a data.frame with User's and Produce's Accuracy or a list for confusion 
#'        matrices.}
#' }
#'
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}},
#' \code{\link[dtwSat]{twdtwApply}}, and 
#' \code{\link[dtwSat]{twdtwClassify}}.
#'
#' @examples
#' \dontrun{
#' load(system.file("lucc_MT/field_samples_ts.RData", package="dtwSat"))
#' set.seed(1)
#' partitions = splitDataset(field_samples_ts, p=0.1, times=5, 
#'                           freq = 8, formula = y ~ s(x, bs="cc"))
#' log_fun = logisticWeight(alpha=-0.1, beta=50) 
#' twdtw_res = lapply(partitions, function(x){
#'    res = twdtwApply(x = x$ts, y = x$patterns, weight.fun = log_fun, n=1)
#'    twdtwClassify(x = res)
#' })
#' cross_validation = twdtwCrossValidation(twdtw_res)
# head(cross_validation, 5)
#' }
NULL

setGeneric("splitDataset", function(object, times, p, ...) standardGeneric("splitDataset"))

#' @rdname Cross-validation
#' @aliases splitDataset
#' @export
setMethod("splitDataset", "twdtwTimeSeries",
          function(object, times=1, p=0.1, ...) splitDataset.twdtwTimeSeries(object, times=times, p=p, ...))

splitDataset.twdtwTimeSeries = function(object, times, p, ...){
  
  partitions = createDataPartition(y = labels(object), times, p, list = TRUE)
  
  res = lapply(partitions, function(I){
    training_ts = subset(object, I)
    validation_ts = subset(object, -I)
    patt = createPatterns(training_ts, ...)
    list(patterns=patt, ts=validation_ts)
  })
  
  res
}

setGeneric("twdtwCrossValidation", function(object, matrix=FALSE) standardGeneric("twdtwCrossValidation"))

#' @rdname Cross-validation
#' @aliases twdtwCrossValidation
#' @export
setMethod("twdtwCrossValidation", "list",
          function(object, matrix) twdtwCrossValidation.twdtwTimeSeries(object, matrix=matrix))

twdtwCrossValidation.twdtwTimeSeries = function(object, matrix){
  
  res = lapply(object, function(x){
    ref = labels(x)$timeseries
    levels = sort(as.character(unique(ref)))
    labels = levels 
    # pred = factor(do.call("rbind", x[])$label, levels, labels)
    pred = do.call("rbind", lapply(x[], function(xx) as.character(xx$label[which.min(xx$distance)])) )
    ref = factor(ref, levels, labels)
    table(Reference=ref, Predicted=pred)
  })
  
  if(!matrix){
    res = do.call("rbind", lapply(seq_along(res), function(i){
      x = res[[i]]
      Users = diag(x) / rowSums(x)
      Producers = diag(x) / colSums(x)
      data.frame(resample=i,label=names(Users), UA = Users, PA = Producers, row.names=NULL)
    }))
  }
  
  res
  
}
