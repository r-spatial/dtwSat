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
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### TWDTW ALIGNMENT


#' @title Perform Time-Weighted Dynamic Time Warping alignment
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the alignments between the temporal patterns and 
#' the time series.
#' 
#' @param ... \link[zoo]{zoo} objects
#' @param patterns a list of \link[zoo]{zoo} objects
#' @param normalize.patterns Normalize queries length. Default is FALSE
#' @param patterns.length An integer. Queries length used with normalize. If not declared 
#' the length of the output queries will be the length of the longest patterns
#' @param x A \link[zoo]{zoo} object with a time series similar 
#' to \code{patterns}. \code{x} must have the same number of attributes
#' and be equal to or longer than the \code{patterns}, 
#' \emph{i.e.} \code{nrow(patterns)<=nrow(x)}
#' @param weight.fun A function to be applied to each element of the local cost matrix 
#' \code{phi} and the time cost matrix \code{psi}. If not declared the time weight is zero, 
#' \emph{i.e.} \code{weight.fun = function(phi, psi) 1*phi + 0*psi}. See 'Details'
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' See \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw} [1]
#' @param n.alignments An integer. The maximun number of alignments to 
#' perform. NULL will return all possible alignments
#' @param span Span between two points of minimum in days, \emph{i.e.} the minimum  
#' interval between two alignments, for details see [2]
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE
#' 
#' @details \code{weight.fun} must have two matrices as arguments. The first is a distance 
#' matrix that comes from the difference in the bands values pont-to-point between the 
#' pattern and the time series. The second matrix is the difference in days pont-to-point 
#' between the pattern and the time series. For a logistic time weight 
#' see \link[dtwSat]{logisticWeight} function, see [3,4] for details. 
#' 
#' @docType methods
#' @return A \code{\link[dtwSat]{twdtw-class}} object
#' 
#' @references 
#' [1] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' @references 
#' [2] M\"uller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited. 
### @references 
### [3] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
### (2015). A Time-Weighted Dynamic Time Warping method for land use and land cover 
### mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
### IEEE Journal of, X, XX-XX.
#' @references 
#' [4] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2015). Open boundary dynamic time warping for satellite image time series classification. 
#' In: IGARSS 2015 2015 IEEE International Geoscience and Remote Sensing Symposium, 2015, Milan. 
#' 2015 IEEE International Geoscience and Remote Sensing Symposium (IGARSS), p. 3349-3352. 
#' 
#' @seealso \code{\link[dtwSat]{twdtw-class}}
#' 
#' @examples
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' 
#' # Perform twdtw analysis for a single pixel 
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'              keep=TRUE)
#' getPatternNames(alig)
#' getAlignments(alig)
#' getMatches(alig)
#' getInternals(alig)
#' 
#' # Perform twdtw for a list pixel 
#' aligs = lapply(template.list, FUN=twdtw, patterns=patterns.list, 
#'                weight.fun = weight.fun, keep=TRUE)
#' 
#' aligs
#' 
#' # Perform twdtw in parallel for a list pixel 
#' # require(parallel)
#' # aligs = mclapply(template.list, FUN=twdtw, patterns=patterns.list, 
#' #                  weight.fun = weight.fun, keep=TRUE, mc.cores=2)
#' # aligs
#'                
#' @export
twdtw =  function(x, ..., patterns=list(...), normalize.patterns=FALSE, 
                  patterns.length=NULL, weight.fun=NULL, dist.method="Euclidean", 
                  step.matrix = symmetric1, n.alignments=NULL, span=0, keep=FALSE)
{

  if(!is(patterns, "list"))
    stop("patterns should be a list of zoo objects")
  if(any(!unlist(lapply(patterns, is.zoo))))
    stop("patterns should be a list of zoo objects")
  if(!is(x, "zoo"))
    stop("x should be of class zoo")
  if(!is(step.matrix, "stepPattern"))
    stop("step.matrix is no stepPattern object")
  if(is.null(weight.fun))
    weight.fun = function(phi, psi) 1*phi 
  if(!is(weight.fun, "function"))
    stop("weight.fun is not a function")
  
  if(normalize.patterns)
    patterns = normalizePatterns(patterns=patterns, patterns.length=patterns.length)
  
  res = .twdtw(x, patterns, weight.fun, dist.method, 
               step.matrix, n.alignments, span, keep)
  res
}



#' @title Apply Time-Weighted Dynamic Time Warping alignment to raster 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function applies a multidimensional Time-Weighted DTW 
#' analysis for each pixel time series. 
#' 
#' @param raster.list A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}.
#' Each layer of the \code{\link[raster]{raster}} object is a time step. See Details
#' @param patterns a list of \link[zoo]{zoo} objects. See \code{\link[dtwSat]{twdtw}} for
#' details
#' @param timeline A vector of class \code{\link[base]{Dates}}. 
#' It must have the length of the layers in the \code{\link[raster]{raster}} object 
#' @param doy A \code{\link[raster]{RasterBrick-class}} or \code{\link[raster]{RasterStack-class}} 
#' with the extent as the objects in \code{raster.list}
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} 
#' object in the format "yyyy-mm-dd"
#' @param by A \code{\link[base]{character}} with the intevals size, \emph{e.g.} ''6 month''
#' @param breaks A vector of class \code{\link[base]{Dates}}
#' @param win.fun A function. A function to be applied to the twdtw results. See Details
#' @param win.size A numeric vector. The size of a processing window in col x row order.
#' Default is a single pixel, \emph{i.e.} \code{win.size=c(1,1)}
#' @param chunk.overlap A numeric vector. The overlap of neighboring chunks of an 
#' raster. Overlap is specified for each dimension of the raster col x row. 
#' Usefull when \code{win.size} bigger than a single pixel. Normally it is the size of 
#' the window. 
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details
#' @param chunk.size An integer. Set the number of cells for each block, 
#' See \code{\link[raster]{blockSize}} for details 
#' @param ... other arguments to pass to the functions \code{\link[dtwSat]{twdtw}} and 
#' \code{win.fun}. Note that the 'win.fun' should explicitly name the arguments
#' 
#' @details ...
#' 
#' @docType methods
#' @return A \code{\link[dtwSat]{twdtw-class}} object
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @examples
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' 
#' # Perform twdtw analysis for a single pixel 
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'              keep=TRUE)
#' getPatternNames(alig)
#' getAlignments(alig)
#' getMatches(alig)
#' getInternals(alig)
#' 
#' # Perform twdtw for a list pixel 
#' #aligs = lapply(X=template, FUN=twdtw, patterns=patterns.list, 
#' #weight.fun = weight.fun, keep=TRUE)
#' 
#' @export
twdtwApply = function(raster.list, patterns, 
                       timeline, doy, breaks=NULL, from=NULL, to=NULL, by=NULL,
                       win.fun = NULL, win.size = c(1,1), chunk.overlap,
                       mc.cores = 1, chunk.size, ...){
  
  # Set time line 
  timeline = as.Date(timeline)
  if(!exists("doy")) {
    array_data = rep(as.numeric(format(timeline, "%j")), each=ncell(raster.list[[1]]))
    doy = array(array_data, dim = dim(raster.list[[1]]))
  }
  raster.list = c(raster.list, doy=doy)
  
  # Set classification intervals 
  if(is.null(breaks))
    breaks = seq(from=as.Date(from), to=as.Date(to), by=by)
  
  # Split and set arguments to functions 
  args = list(..., patterns=patterns, breaks=breaks, from=from, to=to, by=by)
  win_fun = .setFunArgs(fun = win.fun, args = args)
  twdtw_fun = .setFunArgs(fun = twdtw, args = args)
  
  # Set blocks Multi-thread parameters
  minblocks = mc.cores
  if(!missing(chunk.size)) minblocks = round(ncell(raster.list$doy) / chunk.size)
  blocks = blockSize(raster.list$doy, minblocks = minblocks)
  threads = seq(1, blocks$n)
  
  # Set chunk overlap 
  if(missing(chunk.overlap)) chunk.overlap = win.size%/%2
  blocks$r1 = blocks$row - chunk.overlap[1]
  blocks$r1[blocks$r1<1] = 1
  blocks$r2 = blocks$row + blocks$nrows - 1 + chunk.overlap[1]
  blocks$r2[blocks$r2>nrow(raster.list$doy)] = nrow(raster.list$doy)
  
  # Set other arguments 
  bands = names(patterns[[1]])
  names(bands) = bands
  blocks$bands = bands
  blocks$nl = length(breaks) - 1 
  
  fun2 = function(i){
    # Crop block time series 
    array.list = lapply(raster.list, .cropTimeSeries, r1=blocks$r1[i], r2=blocks$r2[i])
    nts = seq(1, ncell(array.list$doy))
    
    # Build zoo time series  
    zoo.list = lapply(nts, .bulidZooFromTSList, x=array.list, 
                      timeline=timeline, bands=blocks$bands)
    
    # Apply TWDTW for each pixel time series 
    twdtw_results = lapply(zoo.list, FUN = twdtw_fun)
    
    # Classify time interval for each pixel 
    res = lapply(twdtw_results, FUN = win_fun)
    
    # aux = lapply(seq(1,120*6,6), function(i) i:(i+5) )
    y = extent(raster.list$doy, r1=blocks$r1[i], r2=blocks$r2[i])
    res_brick = brick(crop(raster.list$doy, y), nl=blocks$nl)
    M = apply(data.frame(res), 1, function(x) as.numeric(x))
    M = array(M, dim = dim(res_brick))
    res_brick = setValues(res_brick, M)
    res_brick
  }
  
  # out.list = lapply(1:2, FUN=fun2)
  out.list = mclapply(threads, FUN=fun2, mc.cores=mc.cores)
  
  out.list$fun = max
  out = do.call(mosaic, out.list)
  out
}


#' @title Logistic weight function 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds a logistic time weight 
#' function to compute the TWDTW local cost matrix
#' 
#' @param alpha numeric. The steepness of logistic weight
#' @param beta numeric. The midpoint of logistic weight
#' @param theta numeric between 0 and 1. The weight of the time 
#' for the TWDTW computation. Use \code{theta=0} to cancel the time weight, 
#' \emph{i.e.} to run the original DTW algorithm. Default is 0.5 
#' 
#' @docType methods
#' @return An \code{\link[base]{expression}} object
#' 
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' weight.fun
#' 
#' @export
logisticWeight = function(alpha, beta, theta=0.5){
  function(phi, psi) (1-theta)*phi + theta*(1 / (1 + exp(1) ^ (alpha * (psi - beta ))))
}


.twdtw =  function(x, patterns, weight.fun, dist.method, 
                   step.matrix, n.alignments, span, keep)
{

  res = lapply(patterns, function(pattern){
    # Adjust columns by name if possible  
    if(!is.null(names(pattern)) & !is.null(names(x)))
      x = x[,names(pattern), drop=FALSE]
    if(ncol(pattern)!=ncol(x))
      stop("Number of columns in pattern and in x don't match.")

    # Get day of the year
    tx = as.numeric(format(index(pattern), "%j"))
    ty = as.numeric(format(index(x), "%j"))

    # Compute local cost matrix 
    phi = dist(pattern, x, method=dist.method)
    # Time cost matrix 
    psi = .g(dist(tx, ty, method=dist.method))
    # Weighted local cost matrix 
    cm = weight.fun(phi, psi)
    
    # Compute cost matris 
    internals = .computecost(cm, step.matrix)
    internals$timeWeight = matrix(psi, nrow = nrow(psi))
    internals$localMatrix = matrix(cm, nrow = nrow(cm))

    # Find low cost candidate paths 
    d = internals$costMatrix[internals$N,1:internals$M]
    endPoints = .findMin(d, index(x), span=span)
    
    # Trace back low cost paths
    if(length(endPoints)<1){
      alignments = NULL
    }else{
      endPoints = endPoints[order(d[endPoints])]
      if(is.null(n.alignments))
        n.alignments = length(endPoints)
      if(length(endPoints) > n.alignments)
        endPoints = endPoints[1:n.alignments]
      # Trace low cost paths (k-th paths)
      matching = .tracepath(dm=internals$directionMatrix, step.matrix=step.matrix, 
                           jmin=endPoints)
      # Get starting point of each path
      startPoints = unlist(lapply(matching, function(map){
        return(map$index2[1])
      }))

      alignments = list()
      alignments$from       = index(x)[startPoints] # This is a vector of Dates
      alignments$to         = index(x)[endPoints] # This is a vector of Dates
      alignments$distance   = d[endPoints] # This is a numeric vector 
      alignments$K          = length(endPoints)
      alignments$matching   = matching # This is a list of data.frames
      if(keep){
        alignments$internals = internals
        alignments$internals$pattern = pattern # Thes is a zoo object
        alignments$internals$x = x # This is a zoo object 
      }
    }
    alignments
  })
  new("twdtw", call=match.call(), alignments=res)
}

.findMin = function(x, timeline, span){
  NonNA = which(!is.na(x))
  dx = diff(x[NonNA])
  index_min = NonNA[which(dx[-length(dx)] < 0 & dx[-1] >= 0)] + 1
  if(tail(dx,1) < 0)
    index_min = c(index_min,length(x))
  order_min = index_min[order(x[index_min])]
  min_out = array()
  for(i in seq_along(index_min)){
    min_out[i] = order_min[i]
    lower_bound = timeline[order_min[i]] - span
    upper_bound = timeline[order_min[i]] + span
    in_span = lower_bound < timeline[order_min] & timeline[order_min] < upper_bound
    order_min[in_span] = NA
  }
  res = min_out[!is.na(min_out)]
  res
}

.g = function(x){
  x[x>(366/2)] = abs(366 - x[x>(366/2)])
  x
}
