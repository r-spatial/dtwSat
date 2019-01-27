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
#' MOD13Q1 250 m 16 days \insertCite{Didan:2015}{dtwSat}. The patterns were built 
#' from ground truth samples of each 
#' crop using Generalized Additive Models (GAM), see \link[dtwSat]{createPatterns}.
#' 
#' @docType data
#' @format A named \code{list} of three \link[zoo]{zoo} objects, ''Soybean'', ''Cotton'', 
#' and ''Maize'', whose indices are \code{\link[base]{Dates}} in the format ''yyyy-mm-dd''.
#' Each node has 6 attributes: ''ndvi'', ''evi'', ''red'', ''nir'', ''blue'', 
#' and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{MOD13Q1.ts},
#' \link[dtwSat]{MOD13Q1.ts.list}, and 
#' \link[dtwSat]{createPatterns}.
#'  
#' @seealso For details about MOD13Q1 see \insertCite{Didan:2015}{dtwSat}.
#' 
#' @references 
#'   \insertAllCited{}
#'   
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}

#' 
"MOD13Q1.patterns.list"


#' @title Data: An example of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a time series based on the 
#' MODIS product MOD13Q1 250 m 16 days \insertCite{Didan:2015}{dtwSat}. 
#' It is an irregularly sampled time series 
#' using the real date of each pixel from ''2009-08-05'' to ''2013-07-31''.
#' 
#' @docType data
#' @format A \link[zoo]{zoo} object, whose indices are \code{\link[base]{Dates}} 
#' in the format ''yyyy-mm-dd''. Each node has 6 attributes: ''ndvi'', 
#' ''evi'', ''red'', ''nir'', ''blue'', and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{MOD13Q1.ts.list},
#' \link[dtwSat]{MOD13Q1.patterns.list}. 
#' 
#' 
#' @seealso For details about MOD13Q1 see \insertCite{Didan:2015}{dtwSat}.
#' 
#' @references 
#'   \insertAllCited{}
#'   
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}

#' 
"MOD13Q1.ts"

#' @title Data: Labels of the satellite time series in MOD13Q1.ts 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description These labels are based on field work.
#' 
#' @docType data
#' @format An object of class \link[base]{data.frame}, whose attributes are: 
#' the label of the crop class ''label'', the start of the crop period ''from'',
#' and the end of the crop period ''to''. The dates are in the format ''yyyy-mm-dd''.
#' 
#' @seealso 
#' \link[dtwSat]{MOD13Q1.ts}. 
#' 
"MOD13Q1.ts.labels"

#' @title Data: A list of satellite time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of time series based on the 
#' MODIS product MOD13Q1 250 m 16 days \insertCite{Didan:2015}{dtwSat}. 
#' It is an irregularly sampled time series 
#' using the real date of each pixel from ''2009-08-05'' to ''2013-07-31''.
#' 
#' @docType data
#' @format A \link[zoo]{zoo} object, whose indices are \code{\link[base]{Dates}} 
#' in the format ''yyyy-mm-dd''. Each node has 6 attributes: ''ndvi'', 
#' ''evi'', ''red'', ''nir'', ''blue'', and ''mir''.
#' 
#' @seealso 
#' \link[dtwSat]{MOD13Q1.ts}, and 
#' \link[dtwSat]{MOD13Q1.patterns.list}. 
#' 
#' 
#' @seealso For details about MOD13Q1 see \insertCite{Didan:2015}{dtwSat}.
#' 
#' @references 
#'   \insertAllCited{}
#'   
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}

#' 
"MOD13Q1.ts.list"

#' @title Data: Pattern time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of patterns with the phenological cycle of: Water,
#' Cotton-Fallow, Forest, Low vegetation, Pasture, Soybean-Cotton, Soybean-Maize, Soybean-Millet, 
#' Soybean-Sunflower, and Wetland. These time series are based on the MODIS product 
#' MOD13Q1 250 m 16 days \insertCite{Didan:2015}{dtwSat}. 
#' The patterns were built from ground truth samples of each 
#' crop using Generalized Additive Models (GAM), see \link[dtwSat]{createPatterns}.
#' 
#' @docType data
#' @format A \link[dtwSat]{twdtwTimeSeries} object.
#' 
#' @seealso For details about MOD13Q1 see \insertCite{Didan:2015}{dtwSat}.
#' 
#' @references 
#'   \insertAllCited{}
#'   
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#' 
"MOD13Q1.MT.yearly.patterns"

