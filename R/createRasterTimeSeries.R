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



#' @title Create raster time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a list of raster time series 
#' 
#' @param x A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}
#' Each layer of the Raster* object is a time step
#' @param timeline A vector of \code{\link[base]{Dates}}
#' It must have the length of the layers in the \code{\link[raster]{RasterBrick-class}} object 
#' @param doy optional. A \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} with the same extent as the 
#' objects in \code{x}
#' @param filepath A character. The path to save the raster time series. If informed the 
#' function saves a raster file for each node in the raster list
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details
#' @param ... other arguments to pass to the functions \code{\link[raster]{writeRaster}} 
#' 
#' @details a list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} objects
#' 
#' @docType methods
#' @return A \code{\link[dtwSat]{twdtw-class}} object
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @examples
#' 
#' ####
#' 
#' @export
createRasterTimeSeries = function(x , timeline, doy, filepath=NULL,
                                  mc.cores = 1, ...){
  
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
