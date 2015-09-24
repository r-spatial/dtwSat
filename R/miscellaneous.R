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
#### UTILITY FUNCTIONS



#' @title Wavelet filter
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a smoothing algorithm to
#' the time series. It computes the a discreat wavelet 
#' smoothing for each dimension in the imput time series.
#' 
#' @param x A \code{zoo} object with the time series
#' @param timeline A vector od dates for the output time series.
#' It must have regular frequency. 
#' @param frequency The frequncy for the output time series
#' @param wf Name of the wavelet filter to use in the decomposition. 
#' Default is "la8"
#' @param J Specifies the depth of the decomposition. This must be a number 
#' less than or equal to log(length(x),2). Default is 1
#' @param boundary Character string specifying the boundary condition. 
#' Default is "periodic". See parameters of \code{\link[waveslim]{mra}}.
#' @param ... see parameters of \code{\link[waveslim]{mra}} in the 
#' packege \pkg{waveslim}
#' @docType methods
#' @return object of class \code{\link[zoo]{zoo}} 
#' 
#' @seealso \link[waveslim]{mra}
#' 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' plot(template$evi, ylab="EVI", xlab="Time")
#' lines(sy$evi, col="red")
#' 
#' ## Plot raw EVI and filtered EVI
#' # require(ggplot2)
#' #df = data.frame(Time=index(template), value=template$evi, variable="Raw")
#' #df = rbind( df, data.frame(Time=index(sy), value=sy$evi, variable="Wavelet filter") )
#' #ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
#' #        geom_line() +
#' #       ylab("EVI")
#'     
#' ## Plot all filter bands
#' # require(ggplot2)
#' # require(reshape2)
#' #df = melt(data.frame(Time=index(sy), sy), id="Time")
#' #ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
#' #       geom_line() 
#'   
#' @export
waveletSmoothing = function(x, timeline=NULL, frequency=NULL, 
                            wf = "la8", J=1, boundary = "periodic", ...)
{
  if( missing(x) )
    stop("Missing either a numeric vector or a zoo object.")
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
  xx = zoo(, timeline)
  xx = merge(x, xx)
  xx = na.approx(xx)
  xx = xx[timeline,]

  # Smoothing 
  df = lapply(as.list(xx), function(d){
    waveslim::mra(x=as.numeric(d), wf=wf, J=J, boundary=boundary)[[paste0("S",J)]]
  })
  res = zoo(data.frame(df), index(xx))
  res
}

