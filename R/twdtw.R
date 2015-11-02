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


#' @title Multidimensional Time-Weighted DTW Alignment
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the alignments of a query within a time series.
#' 
#' @param query A \link[zoo]{zoo} object with a query time series.
#' @param timeseries A \link[zoo]{zoo} object with a time series similar 
#' to \code{query}. \code{timeseries} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' \emph{i.e.} \code{nrow(query)<=nrow(timeseries)}.
#' @param weight A character. ''linear'' for linear weight or ''logistic'' 
#' for logistic weight. Default is NULL that runs the original dtw method,
#' \emph{i.e.} without time weight.
#' @param theta A number. Parameter for ''linear'' time weight. For \code{theta=1} 
#' the time weight is equal to the number of elapsed days.
#' @param alpha A number. The steepness of logistic weight.
#' @param beta A number. The midpoint of logistic weight.
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' See \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw} [2]
#' @param n.alignments An integer. The maximun number of alignments to 
#' perform. NULL will return all possible alignments. 
#' @param span Span between two points of minimum in days, \emph{i.e.} the minimum  
#' interval between two alignments, for details see [1] 
#' @param query.name A query identification
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE
#' @param template is deprecated, please use \code{timeseries} instead
#' @param ...	additional arguments passed to \code{twdtw}
#' 
#' @docType methods
#' @return An object of class \code{\link[dtwSat]{dtwSat-class}}
#' 
#' @references 
#' [1] M\"uller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited. 
#' @references 
#' [2] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' 
#' @seealso \code{\link[dtwSat]{mtwdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' 
#' @examples
#' names(query.list)
#' query.name = "Soybean"
#' alig = twdtw(query.list[[query.name]], timeseries=template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, n.alignments=4, query.name=query.name)
#' alig
#' 
#' @export
twdtw =  function(query, timeseries=NULL, template=NULL,
                  weight=NULL, theta=NULL, alpha=NULL, beta=NULL, 
                  dist.method="Euclidean", step.matrix = symmetric1, 
                  n.alignments=NULL, span=NULL, query.name=NULL, keep=FALSE,
                  ...)
{

  
  if (!missing(template)){
    warning("argument template is deprecated, please use timeseries instead", call. = FALSE)
    timeseries = template
  }
  
  if(!is(query, "zoo"))
    stop("query should be of class zoo")
  if(!is(timeseries, "zoo"))
    stop("timeseries should be of class zoo")
  if(ncol(query)!=ncol(timeseries))
    stop("Number of columns in query and in timeseries don't match")
  if (!is(step.matrix, "stepPattern"))
    stop("step.matrix is no stepPattern object")

  res = .twdtw(query, timeseries, weight, theta, alpha, beta, 
         dist.method, step.matrix, n.alignments, span, 
         query.name, keep, ...)
  res
}

#' @title Performs multiple Time-Weighted DTW 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs the Time-Weighted DTW for a list 
#' of queries
#' 
#' @param query A \link[zoo]{zoo} object with a query time series.
#' @param timeseries A \link[zoo]{zoo} object with a time series similar 
#' to \code{query}. \code{timeseries} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' \emph{i.e.} \code{nrow(query)<=nrow(timeseries)}.
#' @param normalize Normalize queries length. Default is FALSE
#' @param query.length An integer. Queries length used with normalize. If not declared 
#' the length of the output queries will be the length of the longest query
#' @param ... additional arguments passed to \code{\link[dtwSat]{twdtw}}
#' @param template is deprecated, please use \code{timeseries} instead
#' 
#' @docType methods
#' @return An object of class \link[dtwSat]{dtwSat-class} without \code{internals} 
#' or \code{mapping}
#'
#' @seealso \code{\link[dtwSat]{twdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' @examples
#' alig = mtwdtw(query.list, timeseries=template, weight = "logistic", 
#'        alpha = 0.1, beta = 100)
#' alig
#' 
#' @export
mtwdtw = function(query, timeseries=NULL, template=NULL, normalize=FALSE, query.length=NULL, ...){

  
  if (!missing(template)){
    warning("argument template is deprecated, please use timeseries instead", call. = FALSE)
    timeseries = template
  }
  
  if(!is(query, "list"))
    stop("query should be a list of zoo objects")
  
  if(any(!unlist(lapply(query, is.zoo))))
    stop("query should be a list of zoo objects")

  query_names = names(query)
  if(is.null(query_names))
    query_names = seq_along(query)
  
  if(normalize)
    query = normalizeQuery(query=query, query.length=query.length)
  
  res = new("dtwSat")
  res@internals$timeseries = timeseries
  res@alignments = do.call("rbind", lapply(query_names, function(i){
    alig = twdtw(query=query[[i]], timeseries=timeseries, query.name=i, ...)
    getAlignments(alig)
  }))
  res
}


.twdtw =  function(query, timeseries, weight, theta, alpha, beta, 
                   dist.method, step.matrix, n.alignments, span,
                   query.name, keep, ...)
{
  # Align query and timeseries by name if names not null
  if(!is.null(names(query)) & !is.null(names(timeseries)))
    timeseries = timeseries[,names(query)]
  # Local cost
  phi = dist(query, timeseries, method=dist.method)
  # Elapsed time
  psi = matrix(0, nrow = nrow(phi), ncol = ncol(phi))
  if(!is.null(weight)){
    psi = .g(query, timeseries, dist.method)
    psi = switch(weight, 
                 linear   = .linearweight(psi, theta),
                 logistic = .logisticweight(psi, alpha, beta)
    )
  }
  cm = phi + psi
  # Compute cost matris 
  internals = .computecost(cm, step.matrix, ...)
  internals$timeWeight = matrix(psi, nrow = nrow(psi))
  internals$localMatrix = matrix(cm, nrow = nrow(cm))
  internals$query = query
  internals$timeseries = timeseries
  
  d = internals$costMatrix[internals$N,1:internals$M]
  endPoints = .findMin(d, index(timeseries), span=span)
  
  if(length(endPoints)<1){
    alignments = list(quey=numeric(0), from=numeric(0), to=numeric(0), distance=numeric(0))
    mapping = list(index1 = numeric(0), index2 = numeric(0))
  }else{
    endPoints = endPoints[order(d[endPoints])]
    if(is.null(n.alignments))
      n.alignments = length(endPoints)
    if(length(endPoints) > n.alignments)
      endPoints = endPoints[1:n.alignments]
    # Trace low cost paths (k-th paths)
    mapping = .traceback(dm=internals$directionMatrix, step.matrix=step.matrix, 
                         jmin=endPoints, ...)
    # Get starting point of each path
    startPoints = unlist(lapply(mapping, function(map){
      return(map$index2[1])
    }))
    
    if(is.null(query.name))
      query.name = 1
    
    alignments = list(query = query.name,
                      from  = index(timeseries)[startPoints],
                      to    = index(timeseries)[endPoints],
                      distance = d[endPoints],
                      stringsAsFactors = FALSE)
  }
  
  if(keep) return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping, internals=internals))
  
  return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping))
}


.findMin = function(x, timeline, span=NULL){
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

.g = function(query, timeseries, dist.method){
  tx = as.numeric(format(index(query), "%j"))
  ty = as.numeric(format(index(timeseries), "%j"))
  psi = dist(tx, ty, method=dist.method)
  psi[psi>(366/2)] = abs(366 - psi[psi>(366/2)])
  psi
}


.logisticweight = function(x, alpha, beta){
  1 / (1 + exp(1)^(-alpha*(x-beta)))
}

.linearweight = function(x, theta){
  theta * x / (366/2)
}


.traceback = function(dm, step.matrix, jmin, ...){
  n = nrow(dm)
  m = ncol(dm)
  if(is.null(jmin))
    jmin = m
  
  res = lapply(jmin, function(j){
    i = n
    dir = step.matrix
    ns = attr(dir,"npat")
    
    nullrows = dir[,2]==0 & dir[,3]==0
    tmp = dir[!nullrows,,drop=FALSE]
    
    steps.list = lapply(1:ns, function(k){
      sbs = tmp[,1]==k  
      spl = tmp[sbs,-1,drop=FALSE]
      nr = nrow(spl)
      spl[nr:1,,drop=FALSE]
    })
    
    I = c(i)
    J = c(j)
    
    repeat{
      if(i==1)
        break     
      s = dm[i,j]
      if(is.na(s))
        break
      
      steps = steps.list[[s]]
      ns = nrow(steps)
      
      trash = lapply(1:ns, function(k){
        if(i-steps[k,1] > 0){
          I <<- c(i-steps[k,1],I)
          J <<- c(j-steps[k,2],J)
        }   
        NULL
      })
      
      i = I[1]
      j = J[1]
    }
    out = list(index1 = I, index2 = J)
    out
  })
  res
}


.computecost = function(cm, step.matrix, ...){

  cm = rbind(0, cm)
  n = nrow(cm)
  m = ncol(cm)
  lm = matrix(NA, nrow=n, ncol=m)
  lm[1,] = 0
  
  nsteps = dim(step.matrix)[1]
  lm[1,1] = cm[1,1]
  sm = matrix(NA, nrow=n, ncol=m)

  dir = step.matrix
  ns = attr(dir,"npat")
  trash = lapply(1:m, function(j){
    trash = lapply(1:n, function(i){
      if(!is.na(lm[i,j]))
        return(NULL)
      
      clist = numeric(ns)+NA
      for(k in 1:nsteps){
        p = dir[k,1]
        I = i-dir[k,2]
        J = j-dir[k,3]
        if(I>=1 && J>=1) {
          w = dir[k,4]
          if(w == -1) {
            clist[p] = lm[I,J]
          }else{
            clist[p] = clist[p]+w*cm[I,J]
          }
        }
      }
      minc = which.min(clist)
      if(length(minc) > 0){
        lm[i,j] <<- clist[minc]
        sm[i,j] <<- minc
      }
    })
  })
  
  res = list()
  res$costMatrix = lm[-1,]
  res$directionMatrix = sm[-1,]
  res$stepPattern = step.matrix
  res$N = n-1
  res$M = m
  res
}







