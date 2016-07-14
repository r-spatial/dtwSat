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
#   R Package dtwSat - 2016-07-14                             #
#                                                             #
###############################################################

#' @title Time-Weighted Dynamic Time Warping for Satellite Image Time Series 
#' @name dtwSat
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Provides an implementation of the Time-Weighted Dynamic Time Warping 
#' (TWDTW) method for land use and land cover mapping using satellite image time series [1]. 
#' TWDTW is based on the Dynamic Time Warping technique and has achieved high accuracy 
#' for land use and land cover classification using satellite data. The method is based 
#' on comparing unclassified satellite image time series with a set of known temporal 
#' patterns (e.g. phenological cycles associated with the vegetation). Using 'dtwSat' 
#' the user can build temporal patterns for land cover types, apply the TWDTW analysis 
#' for satellite datasets, visualize the results of the time series analysis, produce 
#' land cover maps, and create temporal plots for land cover change analysis.
#' 
#'
#' @references 
#' [1] Maus  V,  Camara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, vol.PP, no.99, pp.1-11.
#' 
#' @seealso \code{\link[dtwSat]{twdtwApply}}
#' 
NULL