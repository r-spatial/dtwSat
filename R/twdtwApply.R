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


#' @title Apply Time-Weighted Dynamic Time Warping to raster time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function applies a multidimensional Time-Weighted DTW 
#' analysis for each pixel location in the raster time series and retrieves 
#' a classified raster time series.  
#' 
#' @param x A list of \code{\link[raster]{Raster-class}}
#' \code{\link[raster]{brick}} or \code{\link[raster]{stack}} object.
#' Each layer of the Raster* object is a time step. 
#' See \code{\link[dtwSat]{buildRasterTimeSeries}}.
#' 
#' @param patterns a list of \link[zoo]{zoo} objects. See \code{\link[dtwSat]{twdtw}} for
#' details.
#' 
#' @param win.fun A function. A function to be applied to the TWDTW results. See Details.
#' 
#' @param win.size A numeric vector. The size of a processing window in col x row order.
#' Default is a single pixel, \emph{i.e.} \code{win.size=c(1,1)}.
#' 
#' @param chunk.overlap A numeric vector. The overlap of neighboring chunks of a 
#' raster. Overlap is specified for each dimension of the raster col x row. 
#' Usefull when \code{win.size} bigger than a single pixel. Normally it is the size of 
#' the window. 
#' 
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details.
#' 
#' @param chunk.size An integer. Set the number of cells for each block, 
#' see \code{\link[raster]{blockSize}} for details.  
#' 
#' @param ... other arguments to pass to the functions \code{\link[dtwSat]{twdtw}} and 
#' \code{win.fun}. Note that the 'win.fun' should explicitly name the arguments.
#' 
#' 
#' @docType methods
#' @return A \code{\link[raster]{brick}} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{buildRasterTimeSeries}}, 
#' \code{\link[dtwSat]{classifyIntervals}}, and  
#' \code{\link[dtwSat]{plotLUCC}}.
#' 
#' 
#' @examples
#' 
#' #####
#'  
#' @export
twdtwApply = function(x, patterns, win.fun, 
                      win.size = c(1,1), 
                      chunk.size = 100, 
                      chunk.overlap,
                      mc.cores = 1, ...){
  
  # Split and set arguments to functions 
  # args = list(weight.fun = weight.fun, breaks = breaks, levels = levels, labels = labels, simplify= TRUE, patterns = patterns.vignette)
  args = list(..., patterns = patterns)
  win_fun = .setFunArgs(fun = win.fun, args = args)
  twdtw_fun = .setFunArgs(fun = twdtw, args = args)
  
  # Set blocks Multi-thread parameters
  minblocks = round(ncell(x$doy) / chunk.size)
  blocks = blockSize(x$doy, minblocks = minblocks)
  threads = seq(1, blocks$n)
  
  # Set chunk overlap 
  if(missing(chunk.overlap)) chunk.overlap = win.size%/%2
  blocks$r1 = blocks$row - chunk.overlap[1]
  blocks$r1[blocks$r1<1] = 1
  blocks$r2 = blocks$row + blocks$nrows - 1 + chunk.overlap[1]
  blocks$r2[blocks$r2>nrow(x$doy)] = nrow(x$doy)
  
  # Reshape raster bands to match pattern bands
  bands = names(patterns[[1]])
  x = x[c(bands,"doy")]
  
  # Get time line 
  timeline = as.Date(names(x$doy), format="date.%Y.%m.%d")
  
  fun2 = function(i){
    # Crop block time series 
    array.list = lapply(x, .cropTimeSeries, r1=blocks$r1[i], r2=blocks$r2[i])
    nts = seq(1, ncell(array.list$doy))
    
    # Build zoo time series  
    zoo.list = lapply(nts, .bulidZoo, x=array.list, timeline=timeline)
    
    # Apply TWDTW for each pixel time series 
    twdtw_results = lapply(zoo.list, FUN = twdtw_fun)
    
    # Classify time interval for each pixel 
    res = lapply(twdtw_results, FUN = win_fun)
    
    # aux = lapply(seq(1,120*6,6), function(i) i:(i+5) )
    M = apply(data.frame(res), 1, function(x) x)
    e = extent(x$doy, r1=blocks$r1[i], r2=blocks$r2[i])
    res_brick = brick(crop(x$doy, e), nl=ncol(M))
    M = array(M, dim = dim(res_brick))
    res_brick = setValues(res_brick, M)
    res_brick
  }
  
  # out.list = lapply(1, FUN=fun2)
  out.list = mclapply(threads, FUN=fun2, mc.cores=mc.cores)
  
  out.list$fun = max
  out = do.call(mosaic, out.list)
  out
}

# Crop raster time series. Returns a 3D array 
.cropTimeSeries = function(x, r1, r2){
  if(is(x, "RasterBrick")){
    y = extent(x, r1, r2)
    x = crop(x, y)
  } else {
    x = x[r1:r2,,]
  }
  alply(x, c(1,2), as.numeric)
}

# Build zoo time series  
.bulidZoo = function(p, x, timeline){
  # Get time series for each band 
  datasets = lapply(x, function(x) x[[p]])
  datasets$doy = getDatesFromDOY(doy=datasets$doy, year=format(timeline, "%Y"))
  idoy = which(names(datasets) %in% c("doy"))
  
  # Remove invalid values 
  k = unlist(lapply(datasets[-idoy], function(x){
    which(x<0|is.na(x))
  }))
  k = c(k, which(duplicated(datasets$doy)))
  if(length(k)>0) datasets = lapply(datasets, function(x) x[-k] )
  
  # Build multi-band zoo object 
  zoo(data.frame(datasets[-idoy]), order.by = datasets$doy)
}

# Match and set a list of arguments to a function 
.setFunArgs = function(fun, ..., args = list(...)){ 
  base_formals = formals(fun)
  base_formals_names = names(base_formals)
  given_formals = args[names(args) %in% base_formals_names]
  missing_formals_names = setdiff(base_formals_names, names(args))
  new_formals = c(base_formals[missing_formals_names], given_formals)
  new_formals = new_formals[base_formals_names]
  formals(fun) = new_formals
  fun
}

