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
#' @param y A \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval
#' @param x A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} objects 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}
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
extractSampleTimeSeries = function(x, y, proj4string = CRS(as.character(NA)), mc.cores = 1){
  
  if(is(y, "data.frame"))
    y = SpatialPointsDataFrame(y[,c("longitude","latitude")], y, proj4string = proj4string)
  
  if(!(is(y, "SpatialPoints") | is(y, "SpatialPointsDataFrame")))
    stop("y is not SpatialPoints")

  # Reproject points to raster projection 
  y = spTransform(y, CRS(projection(x$doy)))

  # res = lapply(s, FUN=.extractTimeSeries, x = x, y = y)
  s = row.names(y)
  names(s) = s
  res = mclapply(s, FUN=.extractTimeSeries, x = x, y = y, mc.cores = mc.cores)
  
  I = which(unlist(lapply(res, is.null)))
  if(length(I)>0)
    warning(paste("the samples ",paste0(s[I], collapse = ","),
                  " are not over the rester extent or they do not overlap the rester time series"), call. = FALSE)
  
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
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param attr A vector character or numeric. The attributes in \code{x} to use 
#' @param freq An integer. The frequency of the output patterns 
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' packege \pkg{mgcv}
#' 
#' 
#' @docType methods
#' 
#' @return A \code{\link[zoo]{zoo}} object 
#' 
#' @examples
#' 
#' ###
#' 
#' 
#' @export
createPattern = function(x, from, to, freq=1, attr, ...){
  
  # Pattern period 
  from = as.Date(from)
  to = as.Date(to)
  
  # Shift dates to match the same period  
  df = do.call("rbind", lapply(x, function(x){
    res = shiftDate(x, year=as.numeric(format(end, "%Y")))
    res = window(res, start = from, end = to)
    data.frame(time=index(res), res)
  }))
  
  dates = as.Date(df$time)
  pred_time = seq(from, to, freq)
  
  fun = function(y, ...){
    df = data.frame(time=as.numeric(dates), y=y)
    fit = gam(data = df, ...)
    predict.gam(fit, newdata = data.frame(time=as.numeric(pred_time)))
  }
  
  if(missing(attr)) attr = names(df)[-1]
  
  res = sapply(as.list(df[attr]), FUN=fun, ...)
  res = zoo(data.frame(res), as.Date(pred_time))
  res
}


.extractTimeSeries = function(p, x, y){
  pto = y[p,]
  # Check if the coordinates are over the raster extent
  I = .getPointsOverRaster(x=x$doy, y=pto)
  if(length(I)<1){
    warning(paste("sample ",p," is not over the rester extent"), call. = FALSE)
    return(NULL)
  }
  # Check if the sample time interval overlaps the raster time series 
  from  = as.Date(pto$from)
  to    = as.Date(pto$to)
  doy   = c(extract(x = x$doy, y = coordinates(pto)))
  year  = format(as.Date(names(x$doy), format="date.%Y.%m.%d"), "%Y")
  timeline = getDatesFromDOY(year = year, doy = doy)
  layer =  which( from - timeline <= 0 )[1]
  nl    =  which(   to - timeline <= 0 )[1] - layer
  if(nl<=0){
    warning(paste("there is no overlap between sample ",p," and the rester time series"), call. = FALSE)
    return(NULL)
  }
  # Extract raster values 
  I = which(names(x) %in% c("doy"))
  ts = data.frame(sapply(x[-I], extract, y=pto, layer = layer, nl = nl))
  zoo(ts, timeline[layer:(layer+nl-1)])
}

.getPointsOverRaster = function(x, y){
  r_extent = extent(x)
  point.x = c(r_extent[2], r_extent[1], r_extent[1], r_extent[2], r_extent[2])
  point.y = c(r_extent[3], r_extent[3], r_extent[4], r_extent[4], r_extent[3])
  r_sp_extent = list(Polygons(list(Polygon(data.frame(x=point.x, y=point.y))), ID="1"))
  r_sp_extent = SpatialPolygons(r_sp_extent, proj4string = CRS(projection(y)))
  res = over(y, r_sp_extent)
  res = row.names(y)[which(!is.na(res))]
  res
}

