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
#   R Package dtwSat - 2016-02-22                             #
#                                                             #
###############################################################

.twdtw = function(x, y, weight.fun, dist.method, step.matrix, 
                  n, span, min.length, theta, keep){
  timeseries=x[[1]]
  res = lapply(as.list(y), function(pattern){
    pattern=pattern[[1]]
    # Adjust columns by name if possible  
    if(!is.null(names(pattern)) & !is.null(names(timeseries)))
      timeseries = timeseries[,names(pattern), drop=FALSE]
    if(ncol(pattern)!=ncol(timeseries))
      stop("Number of columns (attributes) in y and in x don't match.")
    
    # Get day of the year
    ty = index(pattern)
    tx = index(timeseries)
    doyy = as.numeric(format(index(pattern), "%j"))
    doyx = as.numeric(format(index(timeseries), "%j"))
    
    
    # Compute local cost matrix 
    phi = dist(pattern, timeseries, method=dist.method)
    # Time cost matrix 
    psi = .g(dist(doyy, doyx, method=dist.method))
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
      b = .findMin(d, tx, span = span)
      candidates = data.frame(a[b], d[b], b)
    }
    
    # Order maches by similarity 
    I = order(candidates$d)
    if(length(I)<1) return(NULL)
    
    # Sellect alignments 
    if(is.null(n)) n = length(I)
    if(length(I) > n) I = I[1:n]
    
    # Remove overfit 
    I = I[diff(range(ty))*min.length <= tx[candidates$b[I]] - tx[candidates$a[I]]]
    
    alignments = list()
    alignments$from       = tx[candidates$a[I]] # This is a vector of Dates
    alignments$to         = tx[candidates$b[I]] # This is a vector of Dates
    alignments$distance   = candidates$d[I]     # This is a numeric vector 
    alignments$K          = length(I)           # This is an interger 
    
    if(keep){
      # Trace low cost paths (k-th paths)
      matching = .tracepath(dm=internals$directionMatrix, step.matrix=step.matrix, jmin=candidates$b[I])
      alignments$matching = matching    # This is a list of data.frames with the matching points 
      alignments$internals = internals  # These is a list variables used in the TWDTW computation
    }
    alignments
  })
  names(res) = levels(y)
  res
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

#' @useDynLib dtwSat g
.g = function(phi, step.matrix){

  if(!is.loaded("computecost", PACKAGE = "dtwSat", type = "Fortran"))
    stop("Fortran lib is not loaded")

  n = nrow(phi)
  m = ncol(phi)
  res = .Fortran("g", 
      TM = matrix(as.double(phi), n, m),
      N  = as.integer(n),
      M  = as.integer(m),
      PC = as.double(366),
      PACKAGE="dtwSat")
  res$TM
}

