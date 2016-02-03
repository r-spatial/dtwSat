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


###############################################################
#### DATASET DOCUMENTATION


#' @title Data: patterns time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of queries form the crops: Soybean,
#' Cotton, and Maize. These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days. The queries were extracted from more than 100 
#' samples of each class using Generalized Additive Models (GAM), see 
#' package \link[mgcv]{gam}.
#' @docType data
#' @format A named \code{list} of 3 \link[zoo]{zoo} objects, ''Soybean'', ''Cotton'', 
#' and ''Maize'', whose indices are \code{\link[base]{Dates}} in the format ''yyyy-mm-dd''.
#' Each node has 6 attributes: ''ndvi'', ''evi'', ''red'', ''nir'', ''blue'', 
#' and ''mir''.
#' @seealso 
#' \link[dtwSat]{example_ts},
#' \link[dtwSat]{example_ts.list}, and 
#' \link[dtwSat]{createPattern}.
#'  
#' @seealso MOD13Q1 documentation: See 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
#' 
"patterns.list"

#' @title Data: A example of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a time series based on the 
#' MODIS product MOD13Q1 250 m 16 days. It is an irregular time series 
#' using the real date of each pixel from ''2009-08-05'' to ''2013-07-31''.
#' @docType data
#' @format A \link[zoo]{zoo} object, whose indices are \code{\link[base]{Dates}} 
#' in the format ''yyyy-mm-dd''. Each node has 6 attributes: ''ndvi'', 
#' ''evi'', ''red'', ''nir'', ''blue'', and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{example_ts.list},
#' \link[dtwSat]{patterns.list}. 
#' 
#' @seealso The package \pkg{rwtss} provides a client for Web Time 
#' Series Service \url{https://github.com/albhasan/rwtss.git}.
#' 
#' @seealso MOD13Q1 documentation: 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
#' 
"example_ts"



#' @title Data: A list of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of time series based 
#' on the MODIS product MOD13Q1 250 m 16 days. It is an irregular time 
#' series using the real date of each pixel.
#' @docType data
#' @format Each element in the \link[base]{list} has a \link[zoo]{zoo} object, 
#' whose indices are \code{\link[base]{Dates}} in the format ''yyyy-mm-dd''. 
#' Each node has 6 attributes: ''ndvi'', ''evi'', ''red'', ''nir'', ''blue'', 
#' and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{example_ts}, 
#' \link[dtwSat]{patterns.list}. 
#' 
#' @seealso The package \pkg{rwtss} provides a client for Web Time 
#' Series Service \url{https://github.com/albhasan/rwtss.git}.
#' 
#' @seealso MOD13Q1 documentation: 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
#' 
"example_ts.list"


#' @title Data: patterns time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of queries form the crops: ''Cotton-fallow''
#' ''Forest'', ''Soybean-cotton'', ''Soybean-maize'', and ''Soybean-millet''. 
#' These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days. The queries were extracted from field 
#' samples of each class using Generalized Additive Models (GAM), see 
#' package \link[mgcv]{gam}.
#' @docType data
#' @format A named \code{list} of 5 \link[zoo]{zoo} objects, ''Cotton-fallow''
#' ''Forest'', ''Soybean-cotton'', ''Soybean-maize'', and ''Soybean-millet'', 
#' whose indices are \code{\link[base]{Dates}} in the format ''yyyy-mm-dd''.
#' Each node has 6 attributes: ''blue'', ''evi'', ''mir'', ''ndvi'', ''nir'', 
#' and ''red''.
#' @seealso MOD13Q1 documentation: See 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
#' 
"patterns.vignette"



