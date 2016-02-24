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
  I = match(1:3, seq_along(object))
  print(labels(object)[na.omit(I)])
  invisible(NULL)
}

# Show objects of class twdtwMatches 
show.twdtwMatches = function(object){
  cat("An object of class \"twdtwMatches\"\n")
  cat("Number of time series:",length(object@timeseries),"\n")
  cat("Number of Alignments:",length(object),"\n")
  cat("Patterns labels:",as.character(labels(object@patterns)),"\n")
  invisible(NULL)
}

# Show objects of class twdtwRaster 
show.twdtwRaster = function(object){
  cat("An object of class \"twdtwRaster\"\n")
  cat("Time series layers:",names(object),"\n")
  cat("Time range:",paste(min(object@timeline)),"...",paste(max(object@timeline)),"\n")
  cat("dimensions:",dim(object),"(nlayers, nrow, ncol, length)\n")
  cat("resolution:",res(object)," (x, y)\n")
  cat("extent    :",as.vector(extent(object)), "(xmin, xmax, ymin, ymax)\n")
  cat("coord.ref.:",projection(object),"\n") 
  invisible(NULL)
}

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "show", "twdtwTimeSeries",
          definition = show.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod(f = "show", "twdtwMatches",
          definition = show.twdtwMatches)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "show", "twdtwRaster",
          definition = show.twdtwRaster)

