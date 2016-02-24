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
#   R Package dtwSat - 2016-02-18                             #
#                                                             #
###############################################################

#' @title Get time series 
#' @name getTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Generic method to get time series from objects of class twdtw*.
#' 
#' @param object an twdtw* object.
#' @param labels a vector with the time series labels. If not informed the 
#' function retrieves all time series. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @export
setGeneric("getTimeSeries", function(object, ...) standardGeneric("getTimeSeries"))

#' @inheritParams getTimeSeries
#' @describeIn twdtwTimeSeries Get subset of time series from objects of class twdtwTimeSeries.
setMethod("getTimeSeries", "twdtwTimeSeries",
          function(object, labels=NULL) getTimeSeries.twdtwTimeSeries(object=object, labels=labels) )

#' @inheritParams getTimeSeries
#' @describeIn twdtwMatches Get subset of time series from objects of class twdtwMatches.
setMethod("getTimeSeries", "twdtwMatches",
          function(object, labels=NULL) getTimeSeries.twdtwMatches(object=object, labels=labels) )

#' @inheritParams getTimeSeries
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwRaster}}.
#' @param samples a \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval.
#' If missing ''from'' and/or ''to'', their are set to the time range of the object
#' \code{object}. 
#' As additional attribute of \code{samples} can be used as labels for each sample. 
#' See \code{id.labels}. 
#' 
#' @param id.labels a numeric or character with an attribute of \code{samples} to 
#' be used as labels of the samples. Optional.
#' 
#' @param labels character vector with samples labels. It must have one label for each 
#' sample. Optional.
#' 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{samples} is a \code{\link[base]{data.frame}}.
#'
#' @param join.labels a logical. It TRUE the function joins labels that are identical 
#' to a factor. If FALSE a different label is kept for each samples. 
#'
#' @examples 
#' # Creating objects of class twdtwRaster 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, blue, red, nir, mir, timeline=timeline, doy=doy)
#' 
#' # Read field samples 
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", package="dtwSat"), 
#'            what="character")
#' 
#' # Get time series 
#' ts_samples = getTimeSeries(rts, samples=field_samples, proj4string = proj_str, id.labels="label")
#' 
#' @describeIn twdtwRaster Get subsets of time series from objects of class twdtwRaster.
setMethod("getTimeSeries", "twdtwRaster",
          function(object, samples, proj4string = CRS(as.character(NA)), id.labels=NULL, labels=NULL, 
                   join.labels = TRUE){
              if(!is.null(id.labels)) samples$label = as.character(samples[[id.labels]])
              if(!is.null(labels)) samples$label = as.character(labels)
              if(!"label"%in%names(samples)) samples$label = paste0("ts",row.names(samples))
              if(is(samples, "data.frame")){
                if(!is(proj4string, "CRS")) proj4string = try(CRS(proj4string))
                  samples = SpatialPointsDataFrame(samples[,c("longitude","latitude")], samples, proj4string = proj4string)
              }
              if(!(is(samples, "SpatialPoints") | is(samples, "SpatialPointsDataFrame")))
                  stop("samples is not SpatialPoints or SpatialPointsDataFrame")
              extractTimeSeries.twdtwRaster(x=object, y=samples, join.labels=join.labels)
          })
          
extractTimeSeries.twdtwRaster = function(x, y, join.labels){
    
  # Reproject points to raster projection 
  y = spTransform(y, CRS(projection(x)))
  
  # Check if the coordinates are over the raster extent
  pto = .getPointsOverRaster(x, y)
  if(length(pto)<1)
    stop("extents do not overlap")
  if(length(pto)<length(y))
    warning(paste("raster extent does not overlap samples:",paste(pto, collapse = ",")), call. = FALSE)
  
  # Extract time series 
  ts_list = lapply(x, FUN = extract, y = y[pto,])
  
  # Crop period
  if(join.labels){
    res = do.call("joinPatterns", lapply(seq_along(pto), FUN=.extractTimeSeries, pto = pto, x = ts_list, y = y, timeline=index(x)))
  } else {
    res = do.call("joinTimeSeries", lapply(seq_along(pto), FUN=.extractTimeSeries, pto = pto, x = ts_list, y = y, timeline=index(x)))
  }
  res
}

# Get time series from object of class twdtwTimeSeries by labels 
getTimeSeries.twdtwTimeSeries = function(object, labels){
  if(is.null(labels)) labels = labels(object)
  if(is.numeric(labels)) labels = labels(object)[labels]
  I = match(object@labels, labels)
  object@timeseries[na.omit(I)]
}

# Get time series from object of class twdtwMatches by labels 
getTimeSeries.twdtwMatches = function(object, labels) {
  getTimeSeries.twdtwTimeSeries(object@timeseries, labels)
}


.extractTimeSeries = function(pto, p, x, y, timeline){
  s = y[pto[p],]
  # Check if the sample time interval overlaps the raster time series 
  from  = try(as.Date(s$from))
  to    = try(as.Date(s$to))
  label = try(as.character(s$label))
  if(is(label, "try-error")) label = NULL
  doy   = c(x$doy[p,])
  year  = format(timeline, "%Y")
  dates = getDatesFromDOY(year = year, doy = doy)
  if(is.null(from) | is(from, "try-error")) from = dates[1]
  if(is.null(to) | is(to, "try-error")) to = tail(dates,1)
  layer =  which( from - dates <= 0 )[1]
  nl    =  which(   to - dates <= 0 )[1] - layer
  if(nl<=0){
    warning(paste("time period of sample ",pto[p]," does not overlap rester time series"), call. = FALSE)
    return(NULL)
  }
  # Extract raster values 
  I = which(names(x) %in% c("doy"))
  ts = data.frame(sapply(x[-I], function(x) x[p,layer:(layer+nl-1)]))
  dates = dates[layer:(layer+nl-1)]
  k = !duplicated(dates)
  twdtwTimeSeries(timeseries=zoo(ts[k,], dates[k]), labels=label)
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

