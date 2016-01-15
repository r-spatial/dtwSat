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
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################

.onAttach = function(lib, pkg){
  packageStartupMessage(
    sprintf("Loaded dtwSat v%s. See ?dtwSat for help, citation(\"dtwSat\") for use in publication.\n",
            utils::packageDescription("dtwSat")$Version) )
}

#' @import methods
#' @import zoo 
#' @importFrom proxy dist
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @importFrom waveslim mra
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_point geom_path geom_raster geom_polygon geom_bar geom_area theme xlab ylab coord_flip scale_x_continuous scale_y_continuous scale_x_date scale_y_date scale_fill_brewer annotate scale_fill_gradientn scale_fill_manual aes_string waiver facet_wrap ggtitle coord_fixed
#' @importFrom scales pretty_breaks date_format percent
#' @importFrom grDevices terrain.colors gray.colors
#' @importFrom utils tail head 
#' @importFrom plyr alply
#' @importFrom parallel mclapply
#' @importFrom raster blockSize ncell extent crop brick stack setValues mosaic extract projection nlayers writeRaster subset
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPointsDataFrame over CRS spTransform coordinates
#' @importFrom mgcv gam predict.gam
#' @importFrom RColorBrewer brewer.pal
#' 
NULL


### Import and export functions from other packages

#' @importFrom dtw symmetric1
#' @export 
dtw::symmetric1

#' @importFrom dtw symmetric2 
#' @export 
dtw::symmetric2

#' @importFrom dtw asymmetric 
#' @export 
dtw::asymmetric

#' @importFrom dtw rabinerJuangStepPattern 
#' @export 
dtw::rabinerJuangStepPattern
