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
#' see \link[dtwSat]{logisticWeight} function, see [3] for details. 
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
#' @references 
#' [3] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2015). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, X, XX-XX.
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
#' ### Perform twdtw in parallel for a list pixel 
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

