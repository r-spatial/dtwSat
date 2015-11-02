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
#' @import dtw
#' @import zoo
#' @importFrom proxy dist
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @importFrom waveslim mra
#' @importFrom ggplot2 ggplot geom_line geom_point geom_path geom_raster geom_polygon xlab ylab scale_x_continuous scale_y_continuous scale_x_date scale_y_date scale_fill_brewer annotate scale_fill_gradientn aes_string waiver autoplot
#' @importFrom scales pretty_breaks
#' @importFrom grDevices terrain.colors gray.colors
#' @importFrom utils tail head
#'
NULL




