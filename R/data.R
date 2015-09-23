#' @title Query time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a list of queries form the 
#' crops: Soybean, Cotton, and Maize. The time series are based on 
#' the MODIS product MOD13Q1 250 m 16 days. 
#' The queries were exttracted from more than 100 samples of each class
#' using Generalized Additive Models (GAM), see 
#' package \link[mgcv]{gam}.
#' @docType data
#' @format A named list of 3 zoo objects. Each node "Soybean", 
#' "Cotton", and "Maize" has a zoo object with index of class Date
#' and 6 attributes: "ndvi", "evi", "red", "nir", "blue", "mir",
#' and "evi2".
"query.list"


#' @title Template time series
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This dataset has a template time series based on the 
#' MODIS product MOD13Q1 250 m 16 days. See package pkg{rwtss} that is an 
#' R client for Web Time Series Service
#' \url{https://github.com/albhasan/rwtss.git}.
#' @docType data
#' @format A zoo objects with index of class Date and 
#' 6 attributes: "ndvi", "evi", "red", "nir", "blue", "mir",
#' and "evi2".
"template"





