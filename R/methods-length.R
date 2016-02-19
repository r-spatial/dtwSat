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

setGeneric("length")

# Number of time series in objects of class twdtwTimeSeries 
length.twdtwTimeSeries = function(x){
  length(labels(x))
}

# Number of time series in objects of class twdtwMatches
length.twdtwMatches = function(x){
  aligs = getAlignments(x)
  if(length(aligs)<1) return(0)
  nrow(aligs)
}

#' @inheritParams twdtwTimeSeries-class
#' @examples 
#' # Getting the number of time series in objects of class twdtwTimeSeries
#' ex_ts = twdtwTimeSeries(timeseries = example_ts.list)
#' length(ex_ts)
#' 
#' @describeIn twdtwTimeSeries Number of time series in objects of class twdtwTimeSeries
#' @export
setMethod(f = "length", signature = signature("twdtwTimeSeries"), 
          definition = length.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @examples 
#' # Getting the number of time series in objects of class twdtwMatches
#' matches = new("twdtwMatches")
#' length(matches)
#' 
#' @describeIn twdtwMatches Number of time series in objects of class twdtwMatches
#' @export
setMethod(f = "length", signature = signature("twdtwMatches"), 
          definition = length.twdtwMatches)

