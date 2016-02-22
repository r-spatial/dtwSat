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
#   R Package dtwSat - 2016-02-22                             #
#                                                             #
###############################################################

setGeneric("dim")
dim.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = data.frame(as.character(labels(x)), t(sapply(timeseries, dim)))
  names(res) = c("label", "nrow", "ncol")
  res
}

# Number of time series in objects of class twdtwTimeSeries 
setGeneric("length")
length.twdtwTimeSeries = function(x){
  length(labels(x))
}
# Number of time series in objects of class twdtwMatches
length.twdtwMatches = function(x){
  aligs = getAlignments(x)
  if(length(aligs)<1) return(0)
  nrow(aligs)
}

nrow.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = sapply(timeseries, nrow)
  names(res) = as.character(labels(x))
  res
}

ncol.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = sapply(timeseries, ncol)
  names(res) = as.character(labels(x))
  res
}

#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Dimensions of time series from objects of class twdtwTimeSeries.
#' @export 
setMethod(f = "dim", "twdtwTimeSeries",
          definition = dim.twdtwTimeSeries)

#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Number of rows of time series from objects of class twdtwTimeSeries.
#' @export 
setMethod(f = "nrow", "twdtwTimeSeries",
          definition = nrow.twdtwTimeSeries)
          
#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Number of columns of time series from objects of class twdtwTimeSeries.
#' @export 
setMethod(f = "ncol", "twdtwTimeSeries",
          definition = ncol.twdtwTimeSeries)

#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Number of time series in objects of class twdtwTimeSeries.
#' @export
setMethod(f = "length", signature = signature("twdtwTimeSeries"), 
          definition = length.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @describeIn twdtwMatches Number of time series in objects of class twdtwMatches
#' @export
setMethod(f = "length", signature = signature("twdtwMatches"), 
          definition = length.twdtwMatches)

