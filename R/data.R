###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Münster (WWU), Germany                  #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### DATASET DOCUMENTATION


#' @title Data: query time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of queries form the crops: Soybean,
#' Cotton, and Maize. These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days. The queries were exttracted from more than 100 
#' samples of each class using Generalized Additive Models (GAM), see 
#' package \link[mgcv]{gam}.
#' @docType data
#' @format A named \code{list} of 3 \link[zoo]{zoo} objects, "Soybean", "Cotton", 
#' and "Maize", whose indices are Dates in the format "yyyy-mm-dd".
#' Each node has 6 attributes: "ndvi", "evi", "red", "nir", "blue", 
#' "mir", and "evi2".
#' @seealso \link[dtwSat]{template}
#' @seealso MOD13Q1 documentation: See \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
"query.list"


#' @title Data: template time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a template time series based on the 
#' MODIS product MOD13Q1 250 m 16 days. It is an irregular time series 
#' using the real day of each pixel from "2009-08-05" "2013-07-31".
#' @docType data
#' @format A \link[zoo]{zoo} object, whose indices are Dates in the 
#' format "yyyy-mm-dd". Each node has 6 attributes: "ndvi", 
#' "evi", "red", "nir", "blue", "mir", and "evi2".
#' @seealso \link[dtwSat]{query.list}
#' @seealso Package \pkg{rwtss} provides a client for Web Time 
#' Series Service \url{https://github.com/albhasan/rwtss.git}.
#' @seealso MOD13Q1 documentation: \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
"template"





