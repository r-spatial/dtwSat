###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Münster (WWU), Germany                  #
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
#' @param template A \link[zoo]{zoo} object with a template time series similar 
#' to \code{query}. The \code{template} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' \emph{i.e.} \code{nrow(query)<=nrow(template)}.
#' @param weight A character. ''linear'' for linear weight or ''logistic'' 
#' for logistic weight. Default is NULL that runs the original dtw method,
#' \emph{i.e.} without time weight.
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' See \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' @param theta A number. Parameter for ''linear'' time weight. For \code{theta=1} 
#' the time weight is equal to the number of elapsed days.
#' @param alpha A number. The steepness of logistic weight.
#' @param beta A number. The midpoint of logistic weight.
#' @param alignments An integer. The maximun number of alignments to 
#' perform. NULL will return all possible alignments. 
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw}
#' @param window.function see parameter \code{window.type} in \code{\link[dtw]{dtw}} 
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE
#' @param ... other parameters
#' @docType methods
#' @return An object of class \code{\link[dtwSat]{dtwSat-class}} 
#'  
#' @seealso \code{\link[dtwSat]{mtwdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' 
#' @examples
#' names(query.list)
#' alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 50, alignments=4)
#' alig
#' 
#' @export
twdtw =  function(query, template, weight=NULL, dist.method="Euclidean",
                  theta=NULL, alpha=NULL, beta=NULL, alignments=NULL, 
                  step.matrix = symmetric1, window.function = noWindow,
                  keep=FALSE, ...)
{

  if(!is.zoo(query))
    stop("query should be of class zoo.")
  if(!is.zoo(template))
    stop("template should be of class zoo")
  if(ncol(query)!=ncol(template))
    stop("Number of columns in query and in template don't match.")
  if(!is(index(query),"Date"))
    stop("Index in query should be of class Date.")
  if(!is(index(template),"Date"))
    stop("Index in template should be of class Date.")

  .twdtw(query, template, weight, dist.method, theta, alpha, 
         beta, alignments, step.matrix, window.function, keep, ...)
  
}

#' @title Performs multiple Time-Weighted DTW 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description The function performs the Time-Weighted DTW for a list 
#' of queries
#' 
#' @param query A \link[zoo]{zoo} object with a query time series.
#' @param template A \link[zoo]{zoo} object with a template time series similar 
#' to \code{query}. The \code{template} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' @param ... additional arguments passed to \code{\link[dtwSat]{twdtw}}
#' @docType methods
#' @return An object of class \link[base]{data.frame} with alignment attributes
#'  
#' @seealso \code{\link[dtwSat]{twdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' @seealso This function calls a C code to compute a cost matrix. The C 
#' function is internal in the R package \code{\link[dtw]{dtw}} developed 
#' by Toni Giorgino. 
#' @examples
#' alig = mtwdtw(query.list, template, weight = "logistic", 
#'        alpha = 0.1, beta = 50)
#' alig
#' 
#' @export
mtwdtw = function(query, template, ...){
  if(!is.list(query))
    stop("query is not a list.")
  if(!any(unlist(lapply(query, is, "zoo"))))
    stop("query is not a list of zoo objects.")
  
  query.names = names(query)
  if(is.null(query.names))
    query.names = seq_along(query)
  
  res = do.call("rbind", lapply(query.names, function(i){
    data.frame(query=i, twdtw(query[[i]], template)@alignments)
  }))
  
  res
}


.twdtw =  function(query, template, weight=NULL, dist.method="Euclidean",
                  theta=NULL, alpha=NULL, beta=NULL, alignments=NULL, 
                  step.matrix = symmetric1, window.function = noWindow,
                  keep=FALSE, ...)
{

  # Local cost
  delta = proxy::dist(query, template, method=dist.method)
  # Elapsed time
  phi = 0
  if(!is.null(weight)){
    phi = .timeCostMatrix(query, template, dist.method)
    phi = switch(weight, 
                 linear   = .linearweight(phi, theta),
                 logistic = .logisticweight(phi, alpha, beta)
    )
  }
  delta = delta + phi
  
  # Cost matrix
  delta = rbind(0, delta)
  n = nrow(delta)
  m = ncol(delta)
  cm = matrix(NA, nrow=n, ncol=m)
  cm[1,] = 0
  wm = matrix(FALSE, nrow = n, ncol = m)
  wm[window.function(row(wm), col(wm), query.size = n, reference.size = m)] = TRUE
  internals = .computeCM(wm, delta, cm, step.matrix)
  internals$stepPattern = step.matrix
  internals$costMatrix = internals$costMatrix[-1,]
#   internals$directionMatrix = internals$directionMatrix[-1,]
  internals$stepPattern = step.matrix
  internals$N = n-1
  internals$M = m
  internals$query = query
  internals$template = template
  
  # Porform alignments 
  d = internals$costMatrix[internals$N,1:internals$M]
  NonNA = which(!is.na(d))
  diffd = diff(d[NonNA])
  endPoints = NonNA[which(diffd[-length(diffd)] < 0 & diffd[-1] >= 0)] + 1
  if(tail(diffd,1) < 0)
    endPoints = c(endPoints,length(d))
  if( length(endPoints) < 1 ){
    alignments = list(quey=numeric(0),from=numeric(0), to=numeric(0), distance=numeric(0), normalizedDistance=numeric(0))
    mapping = list(index1 = numeric(0), index2 = numeric(0))
  }else{
    endPoints = endPoints[order(d[endPoints])]
    if(is.null(alignments))
      alignments = length(endPoints)
    if(length(endPoints) > alignments)
      endPoints = endPoints[1:alignments]
    # Map low cost paths (k-th paths)
    mapping = lapply(endPoints, function(b){
      return(.kthbacktrack(internals, b))
    })
    
    # Get the starting point of each path
    startPoints = unlist(lapply(mapping, function(map){
      return(map$index2[1])
    }))
    
    alignments = list(from  = startPoints,
                      to    = endPoints,
                      distance           = d[endPoints],
                      normalizedDistance = d[endPoints] / length(query),                      
                      stringsAsFactors = FALSE)
  }
  
  if(keep) return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping, internals=internals))
  
  return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping))
}

.timeCostMatrix = function(query, template, dist.method){ 
  tx = as.numeric(format(index(query), "%j"))
  ty = as.numeric(format(index(template), "%j"))
  phi = proxy::dist(tx, ty, method=dist.method)
  phi[phi>(366/2)] = abs(366 - phi[phi>(366/2)])
  return(phi)
}

.logisticweight = function(x, alpha, beta){
  return( 1 / (1 + exp(1)^(-alpha*(x-beta))) )
}

.linearweight = function(x, theta){
  return( theta * x / 366 )
}


.computeCM = function(...){
  # This function calls a C code to compute the cost mtrix. The C function 
  # is internal in the R package dtw developed by Toni Giorgino. The C code 
  # is included in the src folder of dtwSat Package. 
  .Call("computeCM_Call", PACKAGE="dtwSat", ...)
}


.kthbacktrack = function(alignment, jmin=NULL){
  # This function was adaptaded from the internal function 
  # 'backtrack' in the R package dtw developed by Toni Giorgino
  
  dir = alignment$stepPattern
  npat = attr(dir,"npat")
  
  n = nrow(alignment$costMatrix)
  m = ncol(alignment$costMatrix)
  
  i = n
  j = jmin
  if(is.null(jmin))
    j = alignment$jmin
  
  nullrows = dir[,2]==0 & dir[,3]==0
  tmp = dir[!nullrows,,drop=FALSE]
  
  stepsCache = list()  
  for(k in 1:npat) {
    sbs = tmp[,1]==k  
    spl = tmp[sbs,-1,drop=FALSE]
    nr = nrow(spl)
    stepsCache[[k]] = spl[nr:1,,drop=FALSE]
  }
  
  ii<-c(i)
  jj<-c(j)
  
  repeat {
    if(i==1)
      break	
    s = alignment$directionMatrix[i,j]
    if(is.na(s))
      break
    
    steps = stepsCache[[s]]
    ns = nrow(steps)
    
    for(k in 1:ns) {
      if(i-steps[k,1] > 0) {
        ii = c(i-steps[k,1],ii)
        jj = c(j-steps[k,2],jj)
      }                         
    }
    
    i = ii[1]
    j = jj[1]
  }
  
  internals = list()
  internals$index1 = ii
  internals$index2 = jj
  return(internals)
}







