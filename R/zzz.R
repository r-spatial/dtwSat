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
  
  ## Register TWDTW as a distance function into package proxy
  is_there <- c("TWDTW","twdtw") %in% proxy::pr_DB$get_entry_names()
  sapply(c("TWDTW","twdtw")[is_there], proxy::pr_DB$delete_entry)
  
  proxy::pr_DB$set_entry(FUN   = .twdtwDist,
                         names = c("TWDTW","twdtw"),
                         loop  = FALSE,
                         type  = "metric",
                         description = "Time-Weighted Dynamic Time Warping",
                         reference   = "Maus V, Camara G, Cartaxo R, Sanchez A, Ramos FM, de Queiroz GR (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover mapping. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9 (8), pp. 3729--3739. <doi: 10.1109/JSTARS.2016.2517118>."  )
}

#' @import zoo 
#' @import raster
#' @import ggplot2
#' @import methods
#' @import rgdal 
#' @import snow 
#' @importFrom proxy dist pr_DB
#' @importFrom reshape2 melt
#' @importFrom scales pretty_breaks date_format percent
#' @importFrom grDevices terrain.colors gray.colors
#' @importFrom plyr alply
#' @importFrom sp Polygon Polygons SpatialPoints SpatialPolygons SpatialPointsDataFrame over CRS spTransform coordinates bbox 
#' @importFrom mgcv gam predict.gam 
#' @importFrom RColorBrewer brewer.pal 
#' @importFrom stats xtabs ave window na.omit sd qnorm 
#' @importFrom lubridate month month<- day day<- year year<-
#' @importFrom caret createDataPartition 
#' @importFrom xtable xtable print.xtable
#' @importFrom utils packageDescription flush.console 
#' @useDynLib dtwSat, .registration = TRUE
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
