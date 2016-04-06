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
#' @description This dataset has a list of patterns with the phenological cycle of: Soybean,
#' Cotton, and Maize. These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days [1]. The patterns were build from ground truth samples of each 
#' crop using Generalized Additive Models (GAM), see \link[dtwSat]{createPatterns}.
#' 
#' @docType data
#' @format A named \code{list} of 3 \link[zoo]{zoo} objects, ''Soybean'', ''Cotton'', 
#' and ''Maize'', whose indices are \code{\link[base]{Dates}} in the format ''yyyy-mm-dd''.
#' Each node has 6 attributes: ''ndvi'', ''evi'', ''red'', ''nir'', ''blue'', 
#' and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{example_ts},
#' \link[dtwSat]{example_ts.list}, and 
#' \link[dtwSat]{createPatterns}.
#'  
#' @seealso MOD13Q1 documentation: See 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}.
#' 
#' @references 
#' [1] Friedl MA, Sulla-Menashe D, Tan B, Schneider A, Ramankutty N, Sibley A, Huang X. (2010).
#' MODIS Collection 5 global land cover: Algorithm refinements and characterization of new
#' datasets. Remote Sensing of Environment, 114(1), 168 182.
#' 
"patterns.list"


#' @title Data: An example of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a time series based on the 
#' MODIS product MOD13Q1 250 m 16 days [1]. It is an irregularly sampled time series 
#' using the real date of each pixel from ''2009-08-05'' to ''2013-07-31''.
#' 
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
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}.
#' 
#' @references 
#' [1] Friedl MA, Sulla-Menashe D, Tan B, Schneider A, Ramankutty N, Sibley A, Huang X. (2010).
#' MODIS Collection 5 global land cover: Algorithm refinements and characterization of new
#' datasets. Remote Sensing of Environment, 114(1), 168 182.
#' 
"example_ts"

#' @title Data: Labels of the satellite time series in example_ts 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This labels are based on field work.
#' 
#' @docType data
#' @format An object of class \link[base]{data.frame}, whose attributas are: 
#' the label of the crop class ''label'', the start of the crop period ''from'',
#' and the end of the crop period ''to''. The dates are in the format ''yyyy-mm-dd''.
#' 
#' @seealso 
#' \link[dtwSat]{example_ts}. 
#' 
"example_ts_labels"

#' @title Data: A list of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of time series based on the 
#' MODIS product MOD13Q1 250 m 16 days [1]. It is an irregularly sampled time series 
#' using the real date of each pixel from ''2009-08-05'' to ''2013-07-31''.
#' 
#' @docType data
#' @format A \link[zoo]{zoo} object, whose indices are \code{\link[base]{Dates}} 
#' in the format ''yyyy-mm-dd''. Each node has 6 attributes: ''ndvi'', 
#' ''evi'', ''red'', ''nir'', ''blue'', and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{example_ts}, and 
#' \link[dtwSat]{patterns.list}. 
#' 
#' @seealso The package \pkg{rwtss} provides a client for Web Time 
#' Series Service \url{https://github.com/albhasan/rwtss.git}.
#' 
#' @seealso MOD13Q1 documentation: 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}.
#' 
#' @references 
#' [1] Friedl MA, Sulla-Menashe D, Tan B, Schneider A, Ramankutty N, Sibley A, Huang X. (2010).
#' MODIS Collection 5 global land cover: Algorithm refinements and characterization of new
#' datasets. Remote Sensing of Environment, 114(1), 168 182.
#' 
"example_ts.list"


#' @title Data: patterns time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of patterns with the phenological cycle of: Water,
#' Cotton-Fallow, Forest, Low vegetation, Pasture, Soybean-Cotton, Soybean-Maize, Soybean-Millet, 
#' Soybean-Sunflower, and Wetland. These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days [1]. The patterns were build from ground truth samples of each 
#' crop using Generalized Additive Models (GAM), see \link[dtwSat]{createPatterns}.
#' 
#' @docType data
#' @format A \link[dtwSat]{twdtwTimeSeries} object.
#' 
#' @seealso MOD13Q1 documentation: See 
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}.
#' 
#' @references 
#' [1] Friedl MA, Sulla-Menashe D, Tan B, Schneider A, Ramankutty N, Sibley A, Huang X. (2010).
#' MODIS Collection 5 global land cover: Algorithm refinements and characterization of new
#' datasets. Remote Sensing of Environment, 114(1), 168 182.
#' 
"yearly_patterns_mt"

