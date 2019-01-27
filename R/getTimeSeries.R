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

setGeneric("getTimeSeries", function(object, ...) standardGeneric("getTimeSeries"))
setGeneric("getPatterns", function(object, ...) standardGeneric("getPatterns"))

#' @title Get time series from twdtw* objects
#' @name getTimeSeries
#' @aliases getPatterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Get time series from objects of class twdtw*.
#' 
#' @param object an object of class twdtw*.
#'
#' @param y a \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval.
#' If missing ''from'' and/or ''to'', they are set to the time range of the 
#' \code{object}. 
#' 
#' @param id.labels a numeric or character with an column name from \code{y} to 
#' be used as sample labels. Optional.
#' 
#' @param labels character vector with time series labels. For signature 
#' \code{\link[dtwSat]{twdtwRaster}} this argument can be used to set the 
#' labels for each sample in \code{y}, or it can be combined with \code{id.labels} 
#' to select samples with a specific label.
#' 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}.
#' 
#' @return An object of class \code{\link[dtwSat]{twdtwTimeSeries}}.
#'
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and 
#' \code{\link[dtwSat]{twdtwMatches-class}}
#'
#' @return a list with TWDTW results or an object \code{\link[dtwSat]{twdtwTimeSeries-class}}. 
#'
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples
#' # Getting time series from objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' getTimeSeries(ts, 2)
#' # Getting time series from objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat = twdtwApply(x=ts, y=patt)
#' getTimeSeries(mat, 2)
#' 
#' ## This example creates a twdtwRaster object and extract time series from it. 
#'
#' # Creating objects of class twdtwRaster with evi and ndvi time series 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, timeline=timeline)
#' 
#' # Location and time range 
#' ts_location = data.frame(longitude = -55.96957, latitude = -12.03864, 
#'                          from = "2007-09-01", to = "2013-09-01")
#' prj_string = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#' 
#' # Extract time series 
#' ts = getTimeSeries(rts, y = ts_location, proj4string = prj_string)
#'  
#' autoplot(ts[[1]], facets = NULL) + xlab("Time") + ylab("Value")
#' 
NULL


#' @aliases getTimeSeries-twdtwTimeSeries
#' @inheritParams getTimeSeries
#' @rdname getTimeSeries 
#' @export
setMethod("getTimeSeries", "twdtwTimeSeries",
          function(object, labels=NULL) getTimeSeries.twdtwTimeSeries(object=object, labels=labels) )

#' @aliases getTimeSeries-twdtwMatches
#' @inheritParams getTimeSeries
#' @rdname getTimeSeries 
#' @export
setMethod("getTimeSeries", "twdtwMatches",
          function(object, labels=NULL) getTimeSeries(object=object@timeseries, labels=labels) )

          
#' @aliases getPatterns-twdtwMatches
#' @inheritParams getTimeSeries
#' @rdname getTimeSeries 
#' @export
setMethod("getPatterns", "twdtwMatches",
          function(object, labels=NULL) getTimeSeries(object=object@patterns, labels=labels) )
          
# Get time series from object of class twdtwTimeSeries by labels 
getTimeSeries.twdtwTimeSeries = function(object, labels){
  res = subset(object, labels)
  res@timeseries
}

#' @aliases getTimeSeries-twdtwRaster
#' @inheritParams getTimeSeries
#' @rdname getTimeSeries 
#' @export
setMethod("getTimeSeries", "twdtwRaster",
          function(object, y, labels=NULL, proj4string = NULL, id.labels=NULL){
          
              y = .adjustLabelID(y, labels, id.labels)

              if(!"from"%in%names(y))
                y$from = as.Date(index(object)[1])
              if(!"to"%in%names(y))
                y$to = as.Date(tail(index(object),1))

              y = .toSpatialPointsDataFrame(y, object, proj4string)

              extractTimeSeries.twdtwRaster(object, y)
          })

extractTimeSeries.twdtwRaster = function(x, y){
    
  # Reproject points to raster projection 
  y = spTransform(y, CRS(projection(x)))
  # Check if the coordinates are over the raster extent
  pto = .getPointsOverRaster(x, y)
  if(length(pto)<1)
    stop("Extents do not overlap")
  if(length(pto)<length(y))
    warning(paste("Raster extent does not overlap samples:",paste(pto, collapse = ",")), call. = FALSE)
  
  # Extract time series 
  ts_list = lapply(as.list(x), FUN = extract, y = y[pto,])
  # Crop period
  res = lapply(seq_along(pto), FUN=.extractTimeSeries, pto = pto, x = ts_list, y = y, timeline=index(x))
  ts_null = !sapply(res, is.null)
  res = res[ts_null]
  labels = as.character(y[pto,]$label)[ts_null]
  twdtwTimeSeries(res, labels=labels)
}

.extractTimeSeries = function(p, pto, x, y, timeline){
  s = y[pto[p],]
  # Check if the sample time interval overlaps with the raster time series 
  from  = try(as.Date(s$from))
  to    = try(as.Date(s$to))
  dates = timeline
  if(any(names(x)=="doy")){
    doy   = c(x$doy[p,])
    year  = format(timeline, "%Y")
    dates = getDatesFromDOY(year = year, doy = doy)
  }
  if(is.null(from) | is(from, "try-error")) from = dates[1]
  if(is.null(to) | is(to, "try-error")) to = tail(dates,1)
  # layer =  which( from - dates <= 0 )[1]
  # nl    =  which(   to - dates <= 0 )[1] - layer
  layer =  which.min( abs(from - dates))
  nl    =  which.min( abs(to   - dates)) - layer
  if(nl<=0){
    warning(paste("Time period of sample ",pto[p]," does not overlap rester time series"), call. = FALSE)
    return(NULL)
  }
  # Extract raster values 
  ts = data.frame(lapply(x[!names(x)%in%c("doy")], function(x) x[p,layer:(layer+nl-1)]))
  dates = dates[layer:(layer+nl-1)]
  k = !duplicated(dates)
  zoo(data.frame(ts[k,, drop=FALSE]), dates[k])
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
