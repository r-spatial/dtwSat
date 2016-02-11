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
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @title Extract time series from raster 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function extracts time series from a list of 
#' raster using the spatial location and the time period. 
#' 
#' @param y A \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval.
#' 
#' @param x A list of Raster* objects \code{\link[raster]{Raster-class}}. To build the 
#' raster list use \code{\link[dtwSat]{buildRasterTimeSeries}}. 
#' 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}.
#' 
#' @param mc.cores The number of cores to use, see \code{\link[parallel]{mclapply}} 
#' for details.
#' 
#' @docType methods
#' 
#' @return A list of \code{\link[zoo]{zoo}} objects, one for each location given by \code{y}.
#' 
#' @examples
#' 
#' \dontrun{
#' ##### In the next example we build, extract and plot a multi-band time series. 
#' 
#' #### In this example we build a multi-band MOD13Q1 raster time series. 
#' #### The 'tif' files in 'lucc_MT/tif' have 999 EVI time series 
#' #### from 2007-01-01 to 2013-12-19, that means 160 points with temporal 
#' #### resolution of 16 days. 
#' data_folder = system.file("lucc_MT/data", package = "dtwSat")
#' blue = brick(paste(data_folder,"blue.tif", sep = "/"))
#' red  = brick(paste(data_folder,"red.tif", sep = "/"))
#' nir  = brick(paste(data_folder,"nir.tif", sep = "/"))
#' mir  = brick(paste(data_folder,"mir.tif", sep = "/"))
#' evi  = brick(paste(data_folder,"evi.tif", sep = "/"))
#' ndvi = brick(paste(data_folder,"ndvi.tif", sep = "/"))
#' dates = scan(paste(data_folder,"timeline", sep = "/"), what = "dates")
#' raster_timeseries = 
#'      buildRasterTimeSeries(blue, red, nir, mir, evi, ndvi, timeline = dates)
#'  
#' #### In the next example we extract and plot a multi-band time series from 
#' #### our raster time series built in the example above. 
#' 
#' ## Location and time range 
#' ts_location = data.frame(longitude = -55.96957, latitude = -12.03864, 
#'                          from = "2007-09-01", to = "2013-09-01")
#'  
#' ## Proj string 
#' crs_string = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#' 
#' ## Extract time series 
#' ts = extractTimeSeries(x = raster_timeseries, y = ts_location, 
#'                        proj4string = crs_string)
#'  
#' autoplot(ts[[1]], facets = NULL) + xlab("Time") + ylab("Value")
#' 
#' }
#' 
#' @export
extractTimeSeries = function(x, y, proj4string = CRS(as.character(NA)), mc.cores = 1){
  
  if(is(y, "data.frame")){
    if(!is(proj4string, "CRS")) proj4string = try(CRS(proj4string))
    y = SpatialPointsDataFrame(y[,c("longitude","latitude")], y, proj4string = proj4string)
  }
  
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

