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
#' @aliases twdtwRaster
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
#' names of objects in "\code{...}". 
#' 
#' @param labels a vector of class \code{\link[base]{character}} with 
#' labels of the values in the Raster* objects. This is 
#' useful for categorical Raster* values of land use classes. 
#' 
#' @param levels a vector of class \code{\link[base]{numeric}} with 
#' levels of the values in the Raster* objects. This is 
#' useful for categorical Raster* values of land use classes. 
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
#'  \item{\code{layers}:}{A vector of class \code{\link[base]{character}} 
#'        with the names of the Raster* objects.}
#'  \item{\code{labels}:}{A vector of class \code{\link[base]{factor}} 
#'        with levels and labels of the values in the Raster* objects. This 
#'        is useful for categorical Raster* values of land use classes.}
#' }
#'
#' @seealso   
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{getTimeSeries}},
#' \code{\link[dtwSat]{twdtwMatches-class}}, and 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}
#'
#' @examples 
#' # Creating new object of class twdtwTimeSeries 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = new("twdtwRaster", timeseries = evi, timeline = timeline)
#' 
NULL
twdtwRaster = setClass(
  Class = "twdtwRaster",
  slots = c(timeseries = "list", timeline="Date", layers = "character", labels = "character", levels="numeric"),
  validity = function(object){
    if(!is(object@timeline, "Date")){
      stop("[twdtwTimeSeries: validation] Invalid timeline object, class different from Date.")
    }else{}
    if(any(!(sapply(object@timeseries, is, "RasterBrick") | sapply(object@timeseries, is, "RasterStack") | sapply(object@timeseries, is, "RasterLayer")))){
      stop("[twdtwRaster: validation] Invalid timeseries object, class different from Raster*.")
    }else{}
    if(!is(object@layers, "character")){
      stop("[twdtwTimeSeries: validation] Invalid layers object, class different from character.")
    }else{}
    if( length(object@layers)>0 & length(object@layers)!=length(object@timeseries) ){
      stop("[twdtwTimeSeries: validation] Invalid length, layers and timeseries do not have the same length.")
    }else{}
    if(!is(object@labels, "character")){
      stop("[twdtwTimeSeries: validation] Invalid labels object, class different from character.")
    }else{}
    if(!is(object@levels, "numeric")){
      stop("[twdtwTimeSeries: validation] Invalid levels object, class different from numeric.")
    }else{}
    if( length(object@labels) != length(object@levels) ){
      stop("[twdtwTimeSeries: validation] Invalid length, labels and levels do not have the same length.")
    }else{}
    lapply(object@timeseries, FUN=compareRaster, object@timeseries[[1]], extent=TRUE, rowcol=TRUE, 
                  crs=TRUE, res=TRUE, orig=TRUE, rotation=TRUE, stopiffalse=TRUE)              
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "twdtwRaster",
  definition = 
    function(.Object, timeseries, timeline, doy, layers, labels, levels){
      
      .Object@timeseries = list(doy=brick(), Layer0=brick())
      .Object@timeline = as.Date(0)
      .Object@labels = as.character()
      .Object@levels = numeric()
      if(!missing(timeseries)){
        if(is(timeseries, "RasterBrick") | is(timeseries, "RasterStack") | is(timeseries, "RasterLayer") ) 
           timeseries = list(timeseries)
        if(is.null(names(timeseries))) names(timeseries) = paste0("Layer",seq_along(timeseries))
        .Object@timeseries = c(.Object@timeseries[1], timeseries)
      }
      if(missing(doy)) doy = creat.doy(.Object@timeseries[[2]], timeline)
      .Object@timeseries$doy = doy 
      if(!missing(layers))
        names(.Object@timeseries) = c("doy", layers)
      .Object@layers = names(.Object@timeseries)
      if(!missing(labels))
        .Object@labels = as.character(labels)
      if(missing(levels))
        levels = seq_along(.Object@labels)
      .Object@levels = as.numeric(levels)
      if(!missing(timeline))
        .Object@timeline = as.Date(timeline)
      validObject(.Object)
      names(.Object@timeline) = paste0("date.",format(.Object@timeline,"%Y.%m.%d"))
      .Object@timeseries = lapply(.Object@timeseries, function(x) { names(x)=names(.Object@timeline); x})
      return(.Object)
  }
)


setGeneric(name = "twdtwRaster",  
          def = function(..., timeline, doy=NULL, layers=NULL, labels=NULL, levels=NULL, filepath=NULL) 
            standardGeneric("twdtwRaster")
)


#' @inheritParams twdtwRaster
#' @describeIn twdtwRaster Create object of class twdtwRaster.
#'
#' @examples 
#' \dontrun{
#' # Creating objects of class twdtwRaster 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' ts_evi = twdtwRaster(evi, timeline=timeline)
#' 
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' rts = twdtwRaster(evi, ndvi, blue, red, nir, mir, timeline=timeline, doy=doy)
#' }
#' @export
setMethod(f = "twdtwRaster",  
          definition = function(..., timeline, doy, layers, labels, levels, filepath){
              arg_names = names(list(...))
              not_named = setdiff(as.character(match.call(expand.dots=TRUE)), as.character(match.call(expand.dots=FALSE)))
              if(is.null(arg_names)){ 
                arg_names = not_named
              } else {
                arg_names[arg_names==""] = not_named[arg_names==""]
              }
              x = list(...)
              names(x) = c(arg_names)
              I = which(sapply(x, is, "RasterBrick") | sapply(x, is, "RasterStack") | sapply(x, is, "RasterLayer"))
              if(length(I) < 1)
                stop("there is no Raster* objects in the list of arguments")
              # Split arguments 
              timeseries = x[I]
              dotargs = x[-I]
              creat.twdtwRaster(timeseries=timeseries, timeline=as.Date(timeline), doy=doy,
                                layers=layers, labels=labels, levels=levels, filepath=filepath, dotargs=dotargs)
          })

creat.twdtwRaster = function(timeseries, timeline, doy, layers, labels, levels, filepath, dotargs){
  
  if(is.null(doy)) doy = creat.doy(timeseries[[1]], timeline)
  
  # Check timeline 
  nl = sapply(c(timeseries, doy), nlayers)
  if(any(nl!=length(timeline)))
    stop("raster objects do not have the same length as the timeline")
  
  res = timeseries
  # Save a single file (complete time series) for each raster attribute 
  if (!is.null(filepath)) {
#     print("Saving raster objects. It might take some minutes depending on the objects size and format.")
    dir.create(filepath, showWarnings = FALSE)
    write(as.character(timeline), file = paste(filepath, "timeline", sep="/"))
    aux = c(doy=doy, res)
    res_brick = lapply(names(aux), function(i){
      filename = paste(filepath, i, sep="/")
      dotargs = c(x = aux[[i]], filename = filename, dotargs)
      r = do.call(writeRaster, dotargs)
      r
    })
    names(res_brick) = names(aux)
    res = res_brick[-1]
    doy = res_brick[[1]]
  }
  if(is.null(layers)) layers = names(res)
  new("twdtwRaster", timeseries = res, timeline = timeline, doy=doy, layers = layers, labels = labels, levels=levels)
}

creat.doy = function(x, timeline){
    array_data = rep(as.numeric(format(as.Date(timeline), "%j")), each=ncell(x))
    e = extent(x)
    brick(array(array_data, dim = dim(x)), xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], crs = projection(x))
}

