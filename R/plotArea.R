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


#' @title Plotting accumulated area
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting time series of accumulated area.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwRaster}}.
#' @param time.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the layers to plot. For plot type ''change'' the minimum length 
#' is two.
#' @param time.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the layers. It must have the same 
#' length as time.levels. Default is NULL.
#' @param class.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the levels of the raster values. Default is NULL. 
#' @param class.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the raster values. It must have the same 
#' length as class.levels. Default is NULL.
#' @param class.colors a set of aesthetic values. It must have the same 
#' length as class.levels. Default is NULL. See 
#' \link[ggplot2]{scale_fill_manual} for details.
#' @param perc if TRUE shows the results in percent of area. Otherwise shows the 
#' area in the map units or km2 for no project raster. Default is TRUE.
#'  
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}},
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotMaps}}, 
#' \code{\link[dtwSat]{plotChanges}}, and
#' \code{\link[dtwSat]{plotDistance}}.
#'  
#' @examples
#' \dontrun{
#' 
#' # Create raster time series
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' 
#' # Read fiels samples 
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", 
#'                 package="dtwSat"), what = "character")
#' 
#' # Split samples for training (10%) and validation (90%) using stratified sampling 
#' library(caret) 
#' set.seed(1)
#' I = unlist(createDataPartition(field_samples$label, p = 0.1))
#' training_samples = field_samples[I,]
#' validation_samples = field_samples[-I,]
#' 
#' # Create temporal patterns 
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))
#' 
#' # Run TWDTW analysis for raster time series 
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' r_twdtw = twdtwApply(x=rts, y=temporal_patterns, weight.fun=log_fun, format="GTiff", 
#'                      overwrite=TRUE)
#'                      
#' # Classify raster based on the TWDTW analysis 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' plotArea(r_lucc)
#' 
#' plotArea(r_lucc, perc=FALSE)
#' 
#' }
#' @export
plotArea = function(x, time.levels=NULL, time.labels=NULL, class.levels=NULL, class.labels=NULL, class.colors=NULL, perc=TRUE){
  plot(x, type="area", time.levels=time.levels, time.labels=time.labels, class.levels=class.levels, class.labels=class.labels, class.colors=class.colors, perc=perc)
}

.plotArea = function(x, time.levels, time.labels, class.levels, class.labels, class.colors, perc){

  df.area = do.call("rbind", lapply(time.levels, .getAreaByClass, x, class.levels, class.labels))
  df.area = data.frame(variable = as.numeric(time.labels), df.area, stringsAsFactors = FALSE)
  names(class.colors) = names(df.area)[-1]
  df.area = melt(df.area, "variable", value.name = "Freq", variable.name = "value")
  df.area$Time = as.numeric(df.area$variable)
  df.area$variable = factor(df.area$variable)
  
  if(perc)
    df.area$Freq = df.area$Freq / (sum(df.area$Freq) / length(time.levels))
  
  x.breaks = pretty_breaks()(range(df.area$Time))
  
  gp = ggplot(data=df.area, aes_string(x="Time", y="Freq", fill="value")) +
    geom_area(position = 'stack') + 
    scale_fill_manual(name="Legend", values = class.colors) + 
    scale_x_continuous(expand = c(0.01, 0), breaks = x.breaks) + 
    theme(legend.position = "bottom",
          panel.background = element_blank()) + 
    ylab("Area")
  
  if(perc){
    gp = gp + scale_y_continuous(expand = c(0, 0), labels = percent)
  } else {
    gp = gp + scale_y_continuous(expand = c(0, 0)) 
  }

  gp 
  
}



