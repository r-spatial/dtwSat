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


#' @title Plotting maps
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting time series of maps.
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
#' # Run TWDTW analysis for raster time series 
#' patt = MOD13Q1.MT.yearly.patterns
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
#' time_interval = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), 
#'                     by="12 month")
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' 
#' r_twdtw = twdtwApply(x=rts, y=patt, weight.fun=log_fun, breaks=time_interval, 
#'           filepath="~/test_twdtw", overwrite=TRUE, format="GTiff", mc.cores=3)
#' 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff")
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
  df.map = data.frame(coordinates(x), x[], stringsAsFactors=FALSE)
  df.map = melt(df.map, id.vars = c("x", "y"))
  df.map$value = factor(df.map$value, levels = class.levels, labels = class.labels)
  df.map$variable = time.labels[match(as.character(df.map$variable), time.levels)]
  
  # > df.area
  # variable          value        Freq Time
  # 1      2008  Cotton-fallow 0.292292292 2008
  # 2      2009  Cotton-fallow 0.345345345 2009
  # 3      2010  Cotton-fallow 0.034034034 2010
  
  df.area = do.call("rbind", lapply(time.levels, .getAreaByClass, x, class.levels, class.labels))
  df.area = data.frame(variable = as.numeric(time.labels), df.area, stringsAsFactors = FALSE)
  df.area = melt(df.area, "variable", value.name = "Freq", variable.name = "value")
  df.area$Time = as.numeric(df.area$variable)
  df.area$variable = factor(df.area$variable)
  # df.area = data.frame(prop.table(xtabs(~ variable + value, df.map), 1))
  # df.area$Time = as.numeric(as.character(df.area$variable))
  if(perc){
    df.area$Freq = df.area$Freq / sum(df.area$Freq) / length(time.levels)
  }
    
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



