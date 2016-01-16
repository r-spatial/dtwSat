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
#   R Package dtwSat - 2016-16-01                             #
#                                                             #
###############################################################


#' @title Wavelet filter
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a smoothing algorithm to
#' the time series. It computes a discreat wavelet 
#' smoothing for each dimension in the imput time series.
#' 
#' @param x A \code{\link[zoo]{zoo}} object with the time series
#' @param timeline A vector of dates for the output time series.
#' It must have a regular frequency. 
#' @param frequency The frequncy for the output time series
#' @param wf Name of the wavelet filter used in the decomposition. 
#' Default is "la8"
#' @param J Specifies the depth of the decomposition. This must be a number 
#' less than or equal to log(length(x),2). Default is 1
#' @param boundary Character string specifying the boundary condition. 
#' Default is "periodic". See parameters of \code{\link[waveslim]{mra}}.
#' @param ... other arguments to pass to the function \code{\link[waveslim]{mra}} in the 
#' packege \pkg{waveslim}
#' 
#' @docType methods
#' @return A \code{\link[zoo]{zoo}} object
#' 
#' @seealso \link[waveslim]{mra}
#' 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' 
#' evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
#' if(require(ggplot2)){
#'    ## Plot raw EVI and filtered EVI
#'    gp = autoplot(evi, facets = NULL)
#'    gp
#'     
#'    ## Plot all bands
#'    gp = autoplot(sy, facets = NULL)
#'    gp
#' }
#' 
#' 
#' @export
waveletSmoothing = function(x, timeline=NULL, frequency=NULL, 
                            wf = "la8", J=1, boundary = "periodic", ...)
{
  
  if(!is(x, "zoo"))
    stop("x is not a zoo object")
  
  if(is.null(timeline)){
    if(is.null(frequency)){
      timeline = index(x)
    }else{
      timeline = seq(min(index(x)), max(index(x)), by=frequency)
    }
  }
  
  # Linear interpolation of gaps
  I = which(!is.na(timeline) & !duplicated(timeline))
  timeline = timeline[I]
  y = zoo(order.by = timeline)
  y = merge(x, y)
  y = na.approx(y)
  y = y[timeline,]
  
  # Smoothing 
  df = lapply(as.list(y), function(d){
    mra(x=as.numeric(d), wf=wf, J=J, boundary=boundary, ...)[[paste0("S",J)]]
  })
  res = zoo(data.frame(df), index(y))
  res
}