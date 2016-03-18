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
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) + 
    theme(legend.position = "bottom") + 
    coord_fixed(ratio = 1) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    ggtitle(time.label)
  gp 
  
}


