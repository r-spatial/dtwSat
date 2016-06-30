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


#' @import zoo 
#' @import raster 
#' @import ggplot2
#' @import methods
#' @importFrom proxy dist
#' @importFrom reshape2 melt
#' @importFrom scales pretty_breaks date_format percent
#' @importFrom grDevices terrain.colors gray.colors
#' @importFrom plyr alply
#' @importFrom parallel mclapply
#' @importFrom sp Polygon Polygons SpatialPoints SpatialPolygons SpatialPointsDataFrame over CRS spTransform coordinates bbox 
#' @importFrom mgcv gam predict.gam 
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom stats xtabs ave window na.omit
#' @importFrom lubridate month month<- day day<- year year<-
#' @importFrom caret createDataPartition 
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
