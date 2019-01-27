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
#' passed to "\code{...}". If not provided the layers are set to the 
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
#' If \code{doy} is not provided then at least one Raster* object must be passed through \code{...}.
#'  
#' @param filepath A character. The path to save the raster time series. If provided the 
#' function saves a raster file for each Raster* object in the list, \emph{i.e} one file 
#' for each time series. This way the function retrieves a list of 
#' \code{\link[raster]{RasterBrick-class}}. It is useful when the time series are 
#' originally stored in separated files. See details. 
#'
#' @param object an object of class twdtwRaster.
#'
#' @param x an object of class twdtwRaster.
#'
#' @details The performance of the functions \code{\link[dtwSat]{twdtwApply}} and 
#' \code{\link[dtwSat]{getTimeSeries}} is improved if the Raster* objects are connected 
#' to files with the whole time series for each attribute. 
#'
#' @section Slots :
#' \describe{
#'  \item{\code{timeseries}:}{A list of multi-layer Raster* objects 
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
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples 
#' # Creating a new object of class twdtwTimeSeries 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = new("twdtwRaster", timeseries = evi, timeline = timeline)
#' 
NULL
setClass(
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
      
      .Object@timeseries = list(Layer0=brick())
      .Object@timeline = as.Date(0)
      .Object@labels = as.character()
      .Object@levels = numeric()
      if(!missing(timeseries)){
        if(is(timeseries, "RasterBrick") | is(timeseries, "RasterStack") | is(timeseries, "RasterLayer") ) 
          timeseries = list(timeseries)
        .Object@timeseries = timeseries
        if(is.null(names(.Object@timeseries)))
          names(.Object@timeseries) = paste0("Layer", seq_along(.Object@timeseries)-1)
      } else {
        if(!missing(layers))
          names(.Object@timeseries) = layers
      }
      if(!missing(doy)) 
        .Object@timeseries = c(doy = doy, .Object@timeseries)
      .Object@layers = names(.Object@timeseries)
      if(!missing(labels))
        .Object@labels = as.character(labels)
      if(missing(levels))
        levels = seq_along(.Object@labels)
      .Object@levels = as.numeric(levels)
      if(!missing(timeline))
        .Object@timeline = as.Date(timeline)
      validObject(.Object)
      names(.Object@timeline) = paste0("date.", format(.Object@timeline,"%Y.%m.%d"))
      .Object@timeseries = lapply(.Object@timeseries, function(x) { names(x)=names(.Object@timeline); x})
      return(.Object)
  }
)


setGeneric(name = "twdtwRaster",  
          def = function(..., timeline, doy=NULL, layers=NULL, labels=NULL, levels=NULL, filepath=NULL) 
            standardGeneric("twdtwRaster")
)


#' @inheritParams twdtwRaster
#' @aliases twdtwRaster-create
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
#' rts = twdtwRaster(doy, evi, ndvi, blue, red, nir, mir, timeline = timeline)
#' }
#' @export
setMethod(f = "twdtwRaster",  
          definition = function(..., timeline, doy, layers, labels, levels){
              arg_names = names(list(...))
              not_named = setdiff(as.character(match.call(expand.dots=TRUE)), as.character(match.call(expand.dots=FALSE)))
              if(is.null(arg_names)){ 
                arg_names = not_named
              } else {
                arg_names[arg_names==""] = not_named[arg_names==""]
              }
              x = list(...)
              names(x) = c(arg_names)
              if(missing(doy)){
                if(any(arg_names %in% "doy")){
                  doy = x[[which(arg_names %in% "doy")]]
                  x = x[which(!(arg_names %in% "doy"))]
                }
              }
              I = which(sapply(x, is, "RasterBrick") | sapply(x, is, "RasterStack") | sapply(x, is, "RasterLayer"))
              if(length(I) < 1)
                stop("There are no Raster* objects in the list of arguments")
              # Split arguments 
              timeseries = x[I]
              dotargs = x[-I]
              creat.twdtwRaster(timeseries=timeseries, timeline=as.Date(timeline), doy=doy,
                                layers=layers, labels=labels, levels=levels, dotargs=dotargs)
          })


creat.twdtwRaster = function(timeseries, timeline, doy, layers, labels, levels, dotargs){
  
  # Check timeline 
  nl = sapply(c(timeseries), nlayers)
  if(!is.null(doy))
    nl = c(nlayers(doy), nl)
  if(any(nl!=length(timeline)))
    stop("Raster objects do not have the same length as the timeline")
  
  res = timeseries
  # Save a single file (complete time series) for each raster attribute 
  # if (filepath != "") {
  #   dir.create(filepath, showWarnings = FALSE)
  #   write(as.character(timeline), file = paste(filepath, "timeline", sep="/"))
  #   aux = res
  #   if(!is.null(doy))
  #     aux = c(doy=doy, res)
  #   res_brick = lapply(names(aux), function(i){
  #     filename = paste(filepath, i, sep="/")
  #     dotargs = c(x = aux[[i]], filename = filename, dotargs)
  #     r = do.call(writeRaster, dotargs)
  #     r
  #   })
  #   names(res_brick) = names(aux)
  #   doy = NULL
  #   res = res_brick
  #   if(any(names(res)=="doy")){
  #     res = res_brick[-1]
  #     doy = res_brick[[1]]
  #   }
  # }
  if(is.null(layers)) layers = names(res)
  if(is.null(doy)) 
    return(new("twdtwRaster", timeseries = res, timeline = timeline, layers = layers, labels = labels, levels=levels))
  new("twdtwRaster", timeseries = res, timeline = timeline, doy=doy, layers = layers, labels = labels, levels=levels)
}

.creat.doy = function(x, timeline){
    array_data = rep(as.numeric(format(as.Date(timeline), "%j")), each=ncell(x))
    e = extent(x)
    brick(array(array_data, dim = dim(x)), xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], crs = projection(x))
}


