
#' @title Plotting maps
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting maps and samples.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwAssessment}}.
#' @param samples a character defining the samples to plot 
#' "correct", "incorrect", "all". Default is "all".
#' @param ... other arguments to pass to \code{\link[dtwSat]{twdtwRaster}}
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwAssessment}},
#' \code{\link[dtwSat]{plotAccuracy}}, and
#' \code{\link[dtwSat]{plotAdjustedArea}}.
#'  
#' @example examples/test_twdtw_raster_analysis.R
#' @export
plotMapSamples = function(x, samples="all", ...){
  .plotMapSamples(x, samples, ...)
}

.plotMapSamples = function(x, samples){
  
  x.sp = switch(samples, 
                all       = x@data,
                correct   = x@data[x@data$Predicted == x@data$Reference, ],
                incorrect = x@data[x@data$Predicted != x@data$Reference, ])

  gp = plot(x@map, type="maps") 
  

  df = data.frame(x.sp) 
  df$variable = gp$data[match(df$Period, gp$data$rast.layer),"variable"] 
  df$variable = as.numeric(format(as.Date(df$to), "%Y"))
  
  gp = gp + geom_point(shape = 1, data = df, aes_string(x = "longitude", y = "latitude")) + 
    scale_shape(solid = FALSE) 
  gp 
  
}
