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


#' @title class "twdtwRaster"
#' @name twdtwRaster-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for set of satellite time series.
#' 
#' @param ... objects of class \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}.
#'
#' @param timeline a vector with the dates of the satellite images 
#' in the format of "YYYY-MM-DD".
#'
#' @param layers a vector with the names of the \code{Raster*} objects 
#' passed to "\code{...}". If not informed the layers are set to the 
#' names of the Raster* objects. 
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
#' @param object an object of class twdtwRaster.
#'
#' @param x an object of class twdtwRaster.
#'
#' @details The performance the functions \code{\link[dtwSat]{twdtwApply}} and 
#' \code{\link[dtwSat]{getTimeSeries}} is improved if the Raster* objects are connected 
#' to files with the whole time series for each attribute. 
#'
#' @section Slots :
#' \describe{
#'  \item{\code{timeseries}:}{A list of multi-layers Raster* objects 
#'        with the satellite image time series.}
#'  \item{\code{timeline}:}{A vector of class \code{\link[base]{date}} 
#'        with dates of the satellite images in \code{timeseries}.}
#'  \item{\code{layers}:}{A vector of class \code{\link[base]{factor}} 
#'        with the names of the Raster* objects.}
#' }
#' 
#' @examples 
#' # Creating new object of class twdtwTimeSeries  
#' 
NULL
twdtwRaster = setClass(
  Class = "twdtwRaster",
  slots = c(timeseries = "list", timeline="Date", layers = "factor"),
  validity = function(object){
    if(any(!(sapply(object@timeseries, is, "RasterBrick") | sapply(object@timeseries, is, "RasterStack")))){
      stop("[twdtwRaster: validation] Invalid timeseries object, class different from Raster*.")
    }else{}
    if(!is(object@timeline, "Date")){
      stop("[twdtwTimeSeries: validation] Invalid timeline object, class different from Date.")
    }else{}
    if(!is(object@layers, "factor")){
      stop("[twdtwTimeSeries: validation] Invalid layers object, class different from character.")
    }else{}
    if( length(object@layers)!=0 & length(object@layers)!=length(object@timeseries) ){
      stop("[twdtwTimeSeries: validation] Invalid length, layers and timeseries do not have the same length.")
    }else{}
    lapply(object@timeseries, FUN=compareRaster, object@timeseries[[1]], extent=TRUE, rowcol=TRUE, 
                  crs=TRUE, res=TRUE, orig=TRUE, rotation=TRUE, stopiffalse=TRUE)              
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "twdtwRaster",
  definition = 
    function(.Object, timeseries, doy, timeline, layers){
      .Object@timeseries = list(brick())
      .Object@layers = factor(NULL)
      .Object@timeline = as.Date(0)
      if(!missing(timeseries)){
        if(is(timeseries, "RasterBrick") | is(timeseries, "RasterStack")) 
           timeseries = list(timeseries)
        .Object@timeseries = timeseries
        .Object@layers = factor(paste0("ts",seq_along(.Object@timeseries)))
      }
      if(!missing(layers))
        .Object@layers = factor(layers)
      if(!missing(timeline))
        .Object@timeline = as.Date(timeline)
      validObject(.Object)
      return(.Object)
  }
)

#' @title Create twdtwRaster object 
#' @name twdtwRaster
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Create object of class twdtwRaster.
#' 
#' @inheritParams twdtwRaster-class
#' 
#' @details The performance the functions \code{\link[dtwSat]{twdtwApply}} and 
#' \code{\link[dtwSat]{getTimeSeries}} is improved if the Raster* objects are connected 
#' to files with the whole time series for each attribute. 
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
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), 
#'            what="date")
#' rts=twdtwRaster(evi, ndvi, blue, red, nir, mir, timeline=timeline, doy=doy)
#' 
#' @export
setGeneric(name = "twdtwRaster",  
          def = function(..., timeline, doy=NULL, layers=NULL, filepath=NULL) 
            standardGeneric("twdtwRaster")
)

#' @inheritParams twdtwRaster
#' @describeIn twdtwRaster Create object of class twdtwRaster.
setMethod(f = "twdtwRaster",  
          definition = function(..., timeline, doy, layers, filepath){
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
              timeseries = x[I]
              dotargs = x[-I]
              creat.twdtwRaster(timeseries=timeseries, timeline=as.Date(timeline), 
                                layers=layers, filepath=filepath, dotargs=dotargs)
          })

creat.twdtwRaster = function(timeseries, timeline, layers, filepath, dotargs){
  
  if(is.null(timeseries$doy)) {
    array_data = rep(as.numeric(format(as.Date(timeline), "%j")), each=ncell(timeseries[[1]]))
    e = extent(timeseries[[1]])
    doy = brick(array(array_data, dim = dim(timeseries[[1]])), 
                xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], crs = projection(timeseries[[1]]))
    timeseries$doy=doy
  }
  
  # Check timeline 
  nl = sapply(timeseries, nlayers)
  layer_names = paste0("date.",timeline)
  if(any(nl!=length(layer_names)))
    stop("raster objects do not have the same length as the timeline")
  
  fun = function(x){ names(x) = layer_names; x }
  res = lapply(timeseries, FUN=fun)
  
  # Save a single file (complete time series) for each raster attribute 
  if (!is.null(filepath)) {
    print("Saving raster objects. It might take some minutes depending on the objects size and format.")
    dir.create(filepath, showWarnings = FALSE)
    write(as.character(timeline), file = paste(filepath, "timeline", sep="/"))
    res_brick = lapply(names(res), function(i){
      filename = paste(filepath, i, sep="/")
      dotargs = c(x = res[[i]], filename = filename, dotargs)
      r = do.call(writeRaster, dotargs)
      names(r) = paste0("date.",timeline)
      r
    })
    names(res_brick) = names(res)
    res = res_brick
  }
  if(is.null(layers)) layers = names(res)
  new("twdtwRaster", timeseries = res, timeline = timeline, layers = layers)
}
