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
#   R Package dtwSat - 2015-11-23                             #
#                                                             #
###############################################################


###############################################################
#### EXTRACT AND BUILD TEMPORAL PATTERNS


#' @title Extract time series from raster 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function extracts the time series from a list of 
#' raster given set of spatial location 
#' 
#' @param x A \code{\link[base]{data.frame}} informing: longitude, 
#' latitude, from, to, and classname, for each sample. This can also be 
#' a \code{\link[sp]{SpatialPointsDataFrame}} informing: from, to, and 
#' classname. 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}
#' @param raster.list A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}.
#' @param timeline A vector of class \code{\link[base]{Dates}}. 
#' It must have the length of the layers in the \code{\link[raster]{RasterBrick-class}} objects 
#' @param doy A \code{\link[raster]{RasterBrick-class}} or \code{\link[raster]{RasterStack-class}} 
#' with the same extent as the objects in \code{raster.list}
#' @param year An integer. An base year to set all time series of the class to the same period 
#' @param nsamples An integer. The number of samples for each class. Default extracts 
#' all samples given bx \code{x}
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details
#' 
#' @docType methods
#' 
#' @return A list of \code{\link[base]{data.frame}} objects, one for each class in \code{x}
#' 
#' @examples
#' 
#' ###
#' 
#' 
#' @export
extractSampleTimeSeries = function(x, proj4string=NA, raster.list, timeline, 
                                   doy, year=2000, nsamples=10, mc.cores = 1){
  
  if(is(x, "data.frame"))
    x = SpatialPointsDataFrame(x[,1:2], data=x[,3:5], proj4string=CRS(proj4string))
  
  if(!is(x, "SpatialPointsDataFrame"))
    stop("x is not SpatialPointsDataFrame")
  
  # Get bands 
  bands = names(raster.list)
  names(bands) = bands
  
  # Set time line 
  timeline = as.Date(timeline)
  if(!exists("doy")) {
    array_data = rep(as.numeric(format(timeline, "%j")), each=ncell(raster.list[[1]]))
    doy = array(array_data, dim = dim(raster.list[[1]]))
  }
  raster.list = c(raster.list, doy=doy)
  
  proj4string = CRS(projection(doy))
  
  # Reproject points to raster projection 
  x = spTransform(x, proj4string)
  
  # Get points over the raster
  I = .getPointsOverRaster(x, doy)
  if(length(I)<1)
    stop("there is no points over the raster")
  xx = x[I,]
  
  # Set parameters 
  names(xx) = c("from", "to","classname")
  data = data.frame(xx)
  pattern_names = unique(data$classname)
  names(pattern_names) = pattern_names
  
  fun = function(p){
    I = which(data$classname==p)
    if(!is.null(nsamples) & nsamples<length(I))
      I = sample(I, nsamples)
    if(length(I)<1) return(NULL)
    do.call("rbind", lapply(I, .extractTimeSeries, x=raster.list, 
                            y = xx, timeline = timeline, bands = bands))
  }
  
  # res = lapply(pattern_names[1], FUN=fun)
  res = mclapply(pattern_names, FUN=fun, mc.cores=mc.cores)
  res
}


#' @title Create temporal patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function fits a gam model and retrieves a smoothed 
#' temporal pattern 
#' 
#' @param x A \code{\link[base]{list}} of \code{\link[base]{data.frame}} such as 
#' retrived by \code{\link[dtwSat]{extractSampleTimeSeries}}. 
#' @param freq An integer. The frequency in days of the output patterns 
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' packege \pkg{mgcv}
#' 
#' 
#' @docType methods
#' 
#' @return A list temporal pattern as \code{\link[zoo]{zoo}} objects, 
#' one for each class in \code{x}
#' 
#' @examples
#' 
#' ###
#' 
#' 
#' @export
createPattern = function(x, freq, mc.cores = 1, ...){
  fun = function(x, ...){
    dates = x$time
    time_range = range(dates)
    pred_x = seq(as.Date(time_range[1]), as.Date(time_range[2]), freq)
    xx = x[,-c(1,2)]
    res = sapply(as.list(xx), function(y){
      df = data.frame(x=as.numeric(dates), y=y)
      fit = gam(data = df, ...)
      # fit = gam(data = df, formula = y ~ s(x), family = gaussian(), gamma = 4)
      predict.gam(fit, newdata=data.frame(x=as.numeric(pred_x)))
    })
    zoo(data.frame(res), pred_x)
  }
  # res = lapply(x, FUN=.fitGAM)
  res = mclapply(x, FUN=fun, mc.cores=mc.cores, ...)
  res
}

