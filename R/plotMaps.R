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
#' 
#' }
#' @export
plotMaps = function(x, time.levels=NULL, time.labels=NULL, class.levels=NULL, class.labels=NULL, class.colors=NULL){
  plot(x, type="maps", time.levels=time.levels, time.labels=time.labels, class.levels=class.levels, class.labels=class.labels, class.colors=class.colors)
}

.plotMaps = function(x, time.levels, time.labels, class.levels, class.labels, class.colors){

  df.map = data.frame(coordinates(x), x[], stringsAsFactors=FALSE)
  df.map = melt(df.map, id.vars = c("x", "y"))
  df.map$value = factor(df.map$value, levels = class.levels, labels = class.labels)
  df.map$variable = time.labels[match(as.character(df.map$variable), time.levels)]
  
  gp = ggplot(data=df.map, aes_string(x="x", y="y")) +
    geom_raster(aes_string(fill="value")) + 
    scale_fill_manual(name="Legend", values = class.colors) + 
    facet_wrap(~variable) + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) + 
    theme(legend.position = "bottom") + 
    coord_fixed(ratio = 1) + 
    xlab("Longitude") + 
    ylab("Latitude")
  gp 
  
}
