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



#' @title Build raster time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a list of raster time series.
#' 
#' @param x A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}
#' Each layer of the Raster* object is a time step.
#' 
#' @param timeline A vector of \code{\link[base]{Dates}}
#' It must have the length of the layers in the Raster* object. 
#' 
#' @param doy optional. A \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} with the same spatial and temporal extent 
#' as the Raster* objects in \code{x}.
#' 
#' @param filepath A character. The path to save the raster time series. If informed the 
#' function saves a raster file for each Raster* object in the list, \emph{i.e} one file 
#' for each time series. This way the function retrieves an list of 
#' \code{\link[raster]{RasterBrick-class}} objects. See details. 
#' 
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details.
#' 
#' @param ... other arguments to pass to the function \code{\link[raster]{writeRaster}}.
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
#' ####
#'  #require(raster)
#'  ## Read EVI raster file and timeline 
#'  #r = brick(system.file("lucc_MT/raster_ts/evi.tif",  package = "dtwSat"))
#'  #time = read.csv(system.file("lucc_MT/timeline.csv",  package = "dtwSat"), 
#'  #        as.is=TRUE)
#'          
#'  ## Build raster time series 
#'  #raster_ts = buildRasterTimeSeries(x = r, timeline = dates$date)
#' 
#' @export
buildRasterTimeSeries = function(x , timeline, doy, filepath=NULL,
                                  mc.cores = 1, ...){
  
  if(is(x, "RasterBrick") | is(x, "RasterStack")) x = list(x)
  
  timeline = as.Date(timeline)
  if(missing(doy)) {
    array_data = rep(as.numeric(format(as.Date(timeline), "%j")), each=ncell(x[[1]]))
    e = extent(x[[1]])
    doy = brick(array(array_data, dim = dim(x[[1]])), 
                xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], crs = projection(x[[1]]))
  }
  x = c(x, doy=doy)
  
  # Check timeline 
  nl = lapply(x, nlayers)
  if(any(nl!=length(timeline)))
    stop("raster objects do not have the same length as the timeline")
  
  fun = function(x){
    names(x) = paste0("date.",timeline)
    x
  }
  
  res = lapply(x, FUN=fun)
  
  # Save a single file for each band 
  if (!is.null(filepath)) {
    dir.create(filepath, showWarnings = FALSE)
    res_brick = mclapply(names(res), mc.cores = mc.cores, function(i){
      filename = paste(filepath, i, sep="/")
      dotargs = c(x = res[[i]], filename = filename, list(...))
      r = do.call(writeRaster, dotargs)
      names(r) = paste0("date.",timeline)
      r
    })
    names(res_brick) = names(res)
    res = res_brick
  }
  
  res
}
