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


#' @title Plotting distance maps
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting TWDTW distance maps.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwRaster}}.
#' @param time.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the layers to plot. For plot type ''change'' the minimum length 
#' is two.
#' @param time.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the layers. It must have the same 
#' length as time.levels. Default is NULL.
#' @param layers A \link[base]{character} or \link[base]{numeric}
#' vector with the layers/bands of the raster time series. 
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}},
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotArea}}, 
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
#' plotDistance(r_twdtw)
#' 
#' }
#' @export
plotDistance = function(x, time.levels=1, time.labels=1, layers=NULL){
  plot(x, type="distance", time.levels=time.levels, time.labels=time.labels, layers=layers)
}

.plotDistance = function(x, layers, labels, time.label){
          
  df.map = data.frame(coordinates(x), x[], stringsAsFactors=FALSE)
  df.map = melt(df.map, id.vars = c("x", "y"))
  df.map$variable = labels[match(as.character(df.map$variable), names(x))]
    
  gp = ggplot(data=df.map, aes_string(x="x", y="y")) +
    geom_raster(aes_string(fill="value")) +
    scale_fill_gradient(name="TWDTW distance", low="blue", high="red") + 
    facet_wrap(~variable) + 
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
    theme(legend.position = "bottom") + 
    coord_fixed(ratio = 1) + 
    xlab("") +
    ylab("") + 
    #xlab("Longitude") + 
    #ylab("Latitude") + 
    ggtitle(time.label)
  gp 
  
}


