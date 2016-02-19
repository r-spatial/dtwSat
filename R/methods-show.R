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

# Show objects of class twdtwTimeSeries 
show.twdtwTimeSeries = function(object){
  cat("An object of class \"twdtwTimeSeries\"\n")
  cat("Slot \"timeseries\" length:",length(object),"\n")
  cat("Slot \"labels\": ")
  print(labels(object))
  invisible(NULL)
}


# Show objects of class twdtwMatches 
show.twdtwMatches = function(object){
  cat("An object of class \"twdtwMatches\"\n")
  cat("Number of time series:",length(object@timeseries),"\n")
  cat("Number of Alignments:",length(object),"\n")
  cat("Patterns labels:",labels(object@patterns),"\n")
  invisible(NULL)
}


#' @inheritParams twdtwTimeSeries-class
#' 
#' @examples 
#' # Showing objects of class twdtwTimeSeries 
#' ex_ts = twdtwTimeSeries(timeseries = example_ts.list, labels = c(1,2))
#' show(ex_ts)
#' 
#' @describeIn twdtwTimeSeries Show object of class twdtwTimeSeries  
#' @export 
setMethod(f = "show", "twdtwTimeSeries",
          definition = show.twdtwTimeSeries)


#' @inheritParams twdtwMatches-class
#' 
#' @examples 
#' # Showing objects of class twdtwMatches 
#' matches = new("twdtwMatches")
#' show(matches)
#' 
#' @describeIn twdtwMatches Show object of class twdtwMatches  
#' @export 
setMethod(f = "show", "twdtwMatches",
          definition = show.twdtwMatches)

