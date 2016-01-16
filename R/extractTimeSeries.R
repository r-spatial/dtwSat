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
#' @param x A list of \code{\link[raster]{Raster-class}}
#' \code{\link[raster]{brick}} or \code{\link[raster]{stack}} objects 
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
extractTimeSeries = function(x, y, proj4string = CRS(as.character(NA)), mc.cores = 1){
  
  if(is(y, "data.frame"))
    y = SpatialPointsDataFrame(y[,c("longitude","latitude")], y, proj4string = proj4string)
  
  if(!(is(y, "SpatialPoints") | is(y, "SpatialPointsDataFrame")))
    stop("y is not SpatialPoints or SpatialPointsDataFrame")
  
  # Reproject points to raster projection 
  y = spTransform(y, CRS(projection(x$doy)))
  
  # Check if the coordinates are over the raster extent
  pto = .getPointsOverRaster(x=x$doy, y=y)
  if(length(pto)<1)
    stop("extents do not overlap")
  if(length(pto)<length(y))
    warning(paste("raster extent does not overlap samples:",paste(pto, collapse = ",")), call. = FALSE)
  
  # Extract time series 
  ts_list = mclapply(x, FUN = extract, y = y[pto,], mc.cores = mc.cores)
  
  # Crop period
  res = lapply(seq_along(pto), FUN=.extractTimeSeries, pto = pto, x = ts_list, y = y)
  names(res) = pto
  
  I = which(unlist(lapply(res, is.null)))
  if(length(I)>0)
    warning(paste("rester extent or time period do not overlap samples:",paste0(pto[I], collapse = ",")), call. = FALSE)
  
  res
}

.extractTimeSeries = function(pto, p, x, y){
  s = y[pto[p],]
  # Check if the sample time interval overlaps the raster time series 
  from  = try(as.Date(s$from))
  to    = try(as.Date(s$to))
  doy   = c(x$doy[p,])
  year  = format(as.Date(names(x$doy[p,]), format="date.%Y.%m.%d"), "%Y")
  timeline = getDatesFromDOY(year = year, doy = doy)
  if(is.null(from) | is(from, "try-error")) from = timeline[1]
  if(is.null(to) | is(to, "try-error")) to = tail(timeline,1)
  layer =  which( from - timeline <= 0 )[1]
  nl    =  which(   to - timeline <= 0 )[1] - layer
  if(nl<=0){
    warning(paste("time period of sample ",pto[p]," does not overlap rester time series"), call. = FALSE)
    return(NULL)
  }
  # Extract raster values 
  I = which(names(x) %in% c("doy"))
  ts = data.frame(sapply(x[-I], function(x) x[p,layer:(layer+nl-1)]))
  timeline = timeline[layer:(layer+nl-1)]
  k = !duplicated(timeline)
  zoo(ts[k,], timeline[k])
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

