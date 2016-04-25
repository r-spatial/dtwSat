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
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @title Get dates from year and day of the year
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the date corresponding to the ginven 
#' year and day of the year.
#' 
#' @param year An vector with the years.
#' @param doy An vector with the day of the year. 
#' It must have the same lenght as \code{year}.
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{Dates}} object.
#' 
#' @seealso \link[dtwSat]{shiftDates} 
#' 
#' @examples
#' year = c(2000, 2001)
#' doy = c(366, 365)
#' dates = getDatesFromDOY(year, doy)
#' dates
#'
#' @export
getDatesFromDOY = function(year, doy){
  res = as.Date(paste(as.numeric(year), as.numeric(doy)), format="%Y %j", origin="1970-01-01")
  I = which(diff(res)<0)+1
  res[I] = as.Date(paste0(as.numeric(format(res[I],"%Y"))+1, format(res[I], "-%m-%d")))
  res
}



#' @title Shift dates 
#' @name shiftDates
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function shifts the dates of the time series to a 
#' given base year. 
#' 
#' @param object \code{\link[dtwSat]{twdtwTimeSeries}} objects, 
#' \code{\link[zoo]{zoo}} objects or a list of \code{\link[zoo]{zoo}} objects.
#' 
#' @param year the base year to shit the time series. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}
#'
#' @return An object of the same class as the input \code{object}. 
#'
#' @export
setGeneric("shiftDates", function(object, year=NULL) standardGeneric("shiftDates"))

#' @rdname shiftDates
#' @aliases shiftDates-twdtwTimeSeries
#' @examples
#' patt = twdtwTimeSeries(patterns.list)
#' npatt = shiftDates(patt, year=2005)
#' index(patt)
#' index(npatt)
#' 
#' @export
setMethod("shiftDates", "twdtwTimeSeries",
          function(object, year) 
            do.call("twdtwTimeSeries", lapply(as.list(object), FUN=shiftDates.twdtwTimeSeries, year=year)))

#' @rdname shiftDates
#' @aliases shiftDates-list
#' @export
setMethod("shiftDates", "list",
          function(object, year) 
            shiftDates(twdtwTimeSeries(object), year=year)[])

#' @rdname shiftDates
#' @aliases shiftDates-zoo
#' @export
setMethod("shiftDates", "zoo",
          function(object, year) 
            shiftDates(twdtwTimeSeries(object), year=year)[[1]])

            
shiftDates.twdtwTimeSeries = function(x, year){
  labels = as.character(labels(x))
  x = x[[1]]
  dates = index(x)
  last_date = tail(dates, 1)
  shift_days = as.numeric(last_date - as.Date(paste0(year,format(last_date, "-%m-%d"))))
  d = as.numeric(dates) - shift_days
  new("twdtwTimeSeries", timeseries=zoo(data.frame(x), as.Date(d)), labels=labels)
}


#' @title Classification assessment 
#' @name Assessment
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This functions create data partitions and compute assessment metrics. 
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
#'  \item{\code{twdtwAssess}:}{The function \code{splitDataset} performs the assessment of 
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
#' assessment = twdtwAssess(twdtw_res)
#' head(assessment, 5)
#' }
NULL

setGeneric("splitDataset", function(object, times, p, ...) standardGeneric("splitDataset"))

#' @rdname Assessment
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


setGeneric("twdtwAssess", function(object, matrix=FALSE) standardGeneric("twdtwAssess"))

#' @rdname Assessment
#' @aliases twdtwAssess
#' @export
setMethod("twdtwAssess", "list",
          function(object, matrix) twdtwAssess.twdtwTimeSeries(object, matrix=matrix))

twdtwAssess.twdtwTimeSeries = function(object, matrix){

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

