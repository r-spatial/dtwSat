#' @title Query time series
#' @description This dataset has a list of queries form the 
#' crops: Soybean, Cotton, and Maize. The time series are based on 
#' the MODIS product MOD13Q1 250 m 16 days. 
#' The queries were exttracted from more than 100 samples of each class
#' using Generalized Additive Models (GAM), see \code{\link[gam]{gam}} 
#' in package \pkg{gam}. 
#' @docType data
#' @author Victor Maus, 2015-08-06
"query.list"


#' @title Template time series
#' @description This dataset has a template time series based on the 
#' MODIS product MOD13Q1 250 m 16 days. See package pkg{rwtss} that is an 
#' R client for Web Time Series Service
#' url{https://github.com/albhasan/rwtss.git}.
#' @docType data
#' @author Victor Maus, 2015-08-06
"template"

