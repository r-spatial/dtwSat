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


#' @title Perform Time-Weighted Dynamic Time Warping 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the matches between the temporal patterns and 
#' the time series [1].
#' 
#' @param ... \link[zoo]{zoo} objects.
#' 
#' @param patterns a list of \link[zoo]{zoo} objects.
#' 
#' @param normalize.patterns Normalize patterns length. Default is FALSE.
#' See \link[dtwSat]{normalizePatterns} for details.
#' 
#' @param patterns.length An integer. Patterns length used with \code{patterns.length}. 
#' If not declared the length of the output patterns will be the length of 
#' the longest pattern.
#' 
#' @param x A \link[zoo]{zoo} object with a time series.
#' 
#' @param weight.fun A function that receives a matrix of time differences in days and 
#' returns a matrix of time-weights. If not declared the time-weight is zero. In this 
#' In this case the function runs the standard version of the dynamic time warping. 
#' See details. 
#' 
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' see \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' 
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw} [2].
#' 
#' @param n An integer. The maximun number of matches to perform. 
#' NULL will return all matches.
#' 
#' @param theta numeric between 0 and 1. The weight of the time 
#' for the TWDTW computation. Use \code{theta=0} to cancel the time-weight, 
#' \emph{i.e.} to run the original DTW algorithm. Default is 0.5, meaning that 
#' the time has the same weight as the curve shape in the TWDTW analysis.
#' 
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE. For \code{\link[dtwSat]{plot-method}} methods use \code{keep=TRUE}.
#' 
#' @param span A number. Span between two matches, \emph{i.e.} the minimum  
#' interval between two matches, for details see [3]. If not declared it removes
#' all overlapping matches of the same pattern. To include overlapping matches 
#' of the same pattern use \code{span=0}.
#' 
#' @docType methods
#' @return A \code{\link[dtwSat]{twdtw-class}} object.
#'  
#' @references 
#' [1] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, X, XX-XX.
#' @references 
#' [2] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' @references 
#' [3] M\"uller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{twdtw} through the argument \code{weight.fun}. This will 
#' add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful to analyse time series with phenological cycles
#' of vegetation that are usually bound to seasons. In previous studies by [1] the 
#' logistic weight had better results than the linear for land cover classification. 
#' See [1] for details about the method. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}},
#' \code{\link[dtwSat]{show-method}},  
#' \code{\link[dtwSat]{length-method}}, 
#' \code{\link[dtwSat]{summary-method}},  
#' \code{\link[dtwSat]{plot-method}}, 
#' \code{\link[dtwSat]{getPatternNames}},
#' \code{\link[dtwSat]{getAlignments}}, 
#' \code{\link[dtwSat]{getMatches}},
#' \code{\link[dtwSat]{getInternals}},
#' \code{\link[dtwSat]{linearWeight}}, 
#' \code{\link[dtwSat]{logisticWeight}}, and 
#' \code{\link[dtwSat]{nmatches}}.
#' 
#' @examples
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' 
#' # Perform twdtw analysis for a single pixel 
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun, 
#'              keep=TRUE)
#' getPatternNames(matches)
#' getAlignments(matches)
#' getMatches(matches)
#' getInternals(matches)
#' 
#' # Perform twdtw for a list pixel 
#' matches = lapply(example_ts.list, FUN=twdtw, patterns=patterns.list, 
#'                weight.fun = log_fun, keep=TRUE)
#' 
#' matches
#' 
#' ### Perform twdtw in parallel for a list pixel 
#' # require(parallel)
#' # matches = mclapply(example_ts.list, FUN=twdtw, patterns=patterns.list, 
#' #                  weight.fun = log_fun, keep=TRUE, mc.cores=2)
#' # matches
#'                
#' @export
twdtw =  function(x, ..., patterns=list(...), normalize.patterns=FALSE, 
                  patterns.length=NULL, weight.fun=NULL, dist.method="Euclidean", 
                  step.matrix = symmetric1, n=NULL, span=NULL, theta = 0.5, keep=FALSE)
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
    weight.fun = function(psi) 0 
  if(!is(weight.fun, "function"))
    stop("weight.fun is not a function")
  
  if(normalize.patterns)
    patterns = normalizePatterns(patterns=patterns, patterns.length=patterns.length)
  
  call = match.call()
  
  res = .twdtw(x, patterns, weight.fun, dist.method, 
               step.matrix, n, span, theta, keep, call)
  res
}

.twdtw =  function(x, patterns, weight.fun, dist.method, 
                   step.matrix, n, span, theta, keep, call)
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
    cm = (1-theta)*phi + theta*weight.fun(psi)
    
    # Compute cost matris 
    internals = .computecost(cm, step.matrix)
    internals$timeWeight = matrix(psi, nrow = nrow(psi))
    internals$localMatrix = matrix(cm, nrow = nrow(cm))
    
    # Find low cost candidates 
    d = internals$costMatrix[internals$N,1:internals$M]
    a = internals$startingMatrix[internals$N,1:internals$M]
    if(is.null(span)){
      candidates   = data.frame(a, d)
      candidates   = candidates[ candidates$d==ave(candidates$d, candidates$a, FUN=min), ]
      candidates$b = as.numeric(row.names(candidates))
    }
    else {
      b = .findMin(d, index(x), span = span)
      candidates = data.frame(a[b], d[b], b)
    }
    
    # Order maches by similarity 
    I = order(candidates$d)
    if(length(I)<1) return(NULL)
    
    # Sellect alignments 
    if(is.null(n)) n = length(I)
    if(length(I) > n) I = I[1:n]
    
    alignments = list()
    alignments$from       = index(x)[candidates$a[I]] # This is a vector of Dates
    alignments$to         = index(x)[candidates$b[I]] # This is a vector of Dates
    alignments$distance   = candidates$d[I]           # This is a numeric vector 
    alignments$K          = length(I)                # This is an interger 
    
    if(keep){
      # Trace low cost paths (k-th paths)
      matching = .tracepath(dm=internals$directionMatrix, step.matrix=step.matrix, jmin=candidates$b[I])
      alignments$internals = internals       # These is a list variables used in the TWDTW computation 
      alignments$internals$pattern = pattern # Thes is a zoo object
      alignments$internals$x = x             # This is a zoo object 
      alignments$matching = matching         # This is a list of data.frames with the matching points 
    }
    alignments
  })
  new("twdtw", call=call, alignments=res)
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
    in_span = lower_bound <= timeline[order_min] & timeline[order_min] <= upper_bound
    order_min[in_span] = NA
  }
  res = min_out[!is.na(min_out)]
  res
}

.removeConcurrent = function(I, startPoints, endPoints, d){
  res = !(I | duplicated(startPoints, fromLast = TRUE))
  J = unlist(lapply(unique(startPoints[!res]), function(i){
    J = which(startPoints==i)
    min_j = rep(FALSE, length(J))
    min_j[which.min(d[endPoints[J]])] = TRUE
    min_j
  }))
  res[!res] = J
  res
}

.g = function(x){
  x[x>(366/2)] = abs(366 - x[x>(366/2)])
  x
}

