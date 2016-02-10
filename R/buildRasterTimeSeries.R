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

#' @title Build raster time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a list of raster time series.
#' 
#' @param ... \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} objects with the 
#' same spatial and temporal extents. Each layer of the Raster* object is a time step.
#' \code{...} also accepts other named arguments to pass to the function 
#' \code{\link[raster]{writeRaster}}.
#' 
#' @param timeline A \code{\link[base]{Dates}} or \code{\link[base]{character}} 
#' vector in format "YYYY-MM-DD". It must have the same length as the the layers in 
#' the Raster* object. 
#' 
#' @param doy A \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} with a sequence of days of the year for each pixel. 
#' \code{doy} must have the same spatial and temporal extents as the Raster* objects passed to \code{...}.
#' If \code{doy} is not informed then at least one Raster* object must be passed through \code{...}.
#'  
#' @param filepath A character. The path to save the raster time series. If informed the 
#' function saves a raster file for each Raster* object in the list, \emph{i.e} one file 
#' for each time series. This way the function retrieves an list of 
#' \code{\link[raster]{RasterBrick-class}}. It is useful when the time series are 
#' originally stores in separated files. See details. 
#' 
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details.
#' 
#' @details The performance the functions \code{\link[dtwSat]{twdtwApply}} and 
#' \code{\link[dtwSat]{extractTimeSeries}} is improved if the Raster* objects are connected 
#' to files with the whole time series for each attribute. 
#' 
#' @docType methods
#' 
#' @return details a list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} objects.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwApply}}, and
#' \code{\link[dtwSat]{plotLUCC}}.
#' 
#' @examples
#' 
#' #### In this example we build a multi-band MOD13Q1 raster time series. 
#' #### The 'tif' files in 'lucc_MT/tif' have 999 EVI time series 
#' #### from 2007-01-01 to 2013-12-19, that means 160 points with temporal 
#' #### resolution of 16 days. 
#' blue = brick(system.file("lucc_MT/timeseries/blue.tif",  package = "dtwSat"))
#' red  = brick(system.file("lucc_MT/timeseries/red.tif",  package = "dtwSat"))
#' nir  = brick(system.file("lucc_MT/timeseries/nir.tif",  package = "dtwSat"))
#' mir  = brick(system.file("lucc_MT/timeseries/mir.tif",  package = "dtwSat"))
#' evi  = brick(system.file("lucc_MT/timeseries/evi.tif",  package = "dtwSat"))
#' ndvi = brick(system.file("lucc_MT/timeseries/ndvi.tif",  package = "dtwSat"))
#' dates = scan(system.file("lucc_MT/timeseries/timeline", package = "dtwSat"), what = "dates")
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
#' 
#' #### The next example creates the raster time series and saves one file for each 
#' #### satellite band. This way the function retrieves a rasterBrick object instead 
#' #### of rasterStack that improves the performance of the raster processing. 
#' #### It is useful when the time series are originally stores in separated files. 
#' 
#' # raster_timeseries = buildRasterTimeSeries(blue, red, nir, mir, evi, ndvi,
#' #                      timeline = dates, filepath = "~/test_fun",
#' #                      format="GTiff", overwrite=TRUE)
#'                      
#' @export
buildRasterTimeSeries = function(..., timeline, doy = NULL, filepath = NULL, mc.cores=1){
  
  arg_names = names(list(...))
  not_named = setdiff(as.character(match.call(expand.dots=TRUE)), as.character(match.call(expand.dots=FALSE)))
  if(is.null(arg_names)){ 
    arg_names = not_named
  } else {
    arg_names[arg_names==""] = not_named[arg_names==""]
  }
  x = list(..., doy)
  names(x) = c(arg_names, "doy")
  I = which(sapply(x, is, "RasterBrick") | sapply(x, is, "RasterStack"))
  if(length(I) < 1)
    stop("there is no Raster* objects in the list of arguments")

  # Split arguments 
  raster_objects = x[I]
  dotargs = x[-I]
  
  timeline = as.Date(timeline)
  if(is.null(doy)) {
    array_data = rep(as.numeric(format(as.Date(timeline), "%j")), each=ncell(raster_objects[[1]]))
    e = extent(raster_objects[[1]])
    doy = brick(array(array_data, dim = dim(raster_objects[[1]])), 
                xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], crs = projection(raster_objects[[1]]))
    raster_objects = c(raster_objects, doy=doy)
  }
  
  # Check timeline 
  nl = sapply(raster_objects, nlayers)
  layer_names = paste0("date.",timeline)
  if(any(nl!=length(layer_names)))
    stop("raster objects do not have the same length as the timeline")
  
  
  fun = function(x){ names(x) = layer_names; x }
  res = lapply(raster_objects, FUN=fun)
  
  # Save a single file (complete time series) for each raster attribute 
  if (!is.null(filepath)) {
    dir.create(filepath, showWarnings = FALSE)
    write(as.character(timeline), file = paste(filepath, "timeline", sep="/"))
    res_brick = mclapply(names(res), mc.cores = mc.cores, function(i){
      filename = paste(filepath, i, sep="/")
      dotargs = c(x = res[[i]], filename = filename, dotargs)
      r = do.call(writeRaster, dotargs)
      names(r) = paste0("date.",timeline)
      r
    })
    names(res_brick) = names(res)
    res = res_brick
  }
  
  res
}
