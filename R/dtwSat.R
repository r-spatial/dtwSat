
#' @title Satellite image time series smoothing
#' 
#' @description This function applies a smoothing algorithm to
#' the satellite image time series. It allows to compute the 
#' discreat wavelet smoothing for different levels. This function 
#' only works for constant frequency data, then informa a constant 
#' frequency time series. If it is not the case, either informe a 
#' constant frequency time line or a frequency.
#' 
#' @param x Either a vector of dates or a zoo object with the 
#' time series.
#' @param y A numeric vector.
#' @param timeline A time line vector for the output.
#' @param frequency A numeric with the frequency for output time series.
#' @param method A character vector for the smoothing methods, either 
#' wavelet and level, eg. c(''wavelet'',''1'').
#' @docType methods
#' @examples
#' sy = timeSeriesSmoothing(template$evi, frequency=16, method=c("wavelet",1))
#' plot(template$evi, xlab="Time", ylab="EVI")
#' lines(sy, col="red")
#' @export
timeSeriesSmoothing = function(x, y=NULL, timeline, frequency,
                               method=c("wavelet",1))
{
  if( missing(x) )
    stop("Missing either a numeric vector or a zoo object.")
  
  if(!is.zoo(x))
  {
    if( length(x)!= length(y) )
      stop("Missing numeric vector. y must be a numeric vector the same size as x.")
    I = which(!is.na(x) & !duplicated(x))
    x = x[I]
    y = y[I]
    x = zoo(y, x)
  }
  
  if(missing(timeline)){
    if(missing(frequency)){
      timeline = index(x)
    }else{
      timeline = seq(index(x)[1], index(x)[length(index(x))], by=frequency)
    }
  }
  
  # Linear interpolation of gaps
  I = which(!is.na(timeline) & !duplicated(timeline))
  timeline = timeline[I]
  template = zoo(, timeline)
  template = merge(x, template)
  template = na.approx(template)
  template = template[timeline,]
    
  # Smoothing 
  if(method[1]=="wavelet"){
    sy = mra(as.numeric(template), J=as.numeric(method[2]) ,boundary = "periodic")$S1
    sty = index(template)[seq_along(sy)]
  } else {
    stop("Missing smoothing method. Please choose between discrete wavelet.")
  }
  
  return(zoo(sy, sty))

}


#' @title Multidimensional Time-Weighted DTW analysis
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves one or more possible alignments of a query within 
#' a time series.
#' 
#' @param query A zoo object with the multidimensional time series.
#' @param template A zoo object with the template time series. 
#' It must be iguel or be equal or longer than the length of the query and 
#' the same number of dimensions.
#' @param weight A character. ''linear'' for linear weight or ''logistic'' 
#' for logistic weight. Default is NULL that runs the original dtw method.
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' See \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}
#' @param theta A real. Parameter for linear weight. It is the slope 
#' of the linear TWDTW. For theta equal 1 the time weight is equal to the 
#' number of elapsed days. Default is NULL.
#' @param alpha A real. The steepness of logistic method. Default is NULL.
#' @param beta A real. The midpoint of logistic method. Default is NULL.
#' @param alignments An integer. The maximun number of alignments to 
#' perform. Default is NULL to return all possible alignment. 
#' @param step.matrix see \code{\link{stepPattern}} in package \pkg{dtw}
#' @param window.function see \code{window.type} in package \pkg{dtw}
#' @docType methods
#' @examples
#' alig = mtwdtw(query.list, template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4)
#' alig
#' @export
twdtw =  function(query, template, weight=NULL, dist.method="Euclidean",
                  theta=NULL, alpha=NULL, beta=NULL, alignments=NULL, 
                  step.matrix = symmetric1, window.function = noWindow)
{

  if(!is.zoo(query))
    stop("Missing zoo object. The query must be a zoo object.")
  if(!is.zoo(template))
    stop("Missing zoo object. The template must be a zoo object.")
  if(ncol(query)!=ncol(template))
    stop("Template must have the same number of columns than the query.")
  
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
  out = .Call("computeCM_Call", PACKAGE="dtw", wm, delta, cm, step.matrix)
  out$call = match.call()
  out$stepPattern = step.matrix
  out$costMatrix = out$costMatrix[-1,]
  out$directionMatrix = out$directionMatrix[-1,]
  out$stepPattern = step.matrix
  out$N = n-1
  out$M = m
  out$query = query
  out$template = template
    
  # Porform alignments 
  d = out$costMatrix[out$N,1:out$M]
  NonNA = which(!is.na(d))
  diffd = diff(d[NonNA])
  endPoints = NonNA[which(diffd[-length(diffd)] < 0 & diffd[-1] >= 0)] + 1
  if(tail(diffd,1) < 0)
    endPoints = c(endPoints,length(d))
  if( length(endPoints) < 1 ){
    out$alignments = NULL
  }else{
    endPoints = endPoints[order(d[endPoints])]
    if(is.null(alignments))
      alignments = length(endPoints)
    if(length(endPoints) > alignments)
      endPoints = endPoints[1:alignments]
    # Map low cost paths (k-th paths)
    mapping = lapply(endPoints, function(b){
      return(.kthbacktrack(out, b))
    })
    
    # Get the starting point of each path
    startPoints = unlist(lapply(mapping, function(map){
      return(map$index2[1])
    }))
    
    # Return the alignments
    out$alignments = data.frame(from  = startPoints,
                                to    = endPoints,
                                distance           = d[endPoints],
                                normalizedDistance = d[endPoints] / length(query),                      
                                stringsAsFactors = FALSE)
    out$mapping = mapping
  }
  
  class(out) = "twdtw"
  return(out)
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


#' @title Plotting Time-Weighted DTW 
#' 
#' @description Methods for plotting Time-Weighted DTW objects 
#' returned by twdtw.
#' 
#' @param x \code{\link[dtwSat]{twdtw}} object
#' @param xlab A character, x axis label
#' @param ylab A character, y axis label
#' @param show.dist Display dtw distance for each alignment 
#' @param ... additional arguments passed to plotting functions 
#' @docType methods
#' @examples
#' alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)
#' plot(alig)
#' @export
plot.twdtw = function(x, ylab="Pattern", xlab="Time series", show.dist=FALSE, ...){
  tx = index(x$template)
  ty = index(x$query)
  m = t(x$costMatrix)
  image(m, col=terrain.colors(100), x=tx, y=ty, xlab=xlab, ylab=ylab, ...)
  contour(m, x=tx, y=ty, add=TRUE)
  for(i in 1:length(x$mapping)){
    lines(tx[x$mapping[[i]]$index2], ty[x$mapping[[i]]$index1] ,col="red", lwd=2)
    if(show.dist)
      text(tail(tx[x$mapping[[i]]$index2],1) , tail(ty[x$mapping[[i]]$index1], 1), round(x$alignments$distance[[i]],2))
  }
}


#' @title Method for Time-Weighted DTW 
#' 
#' @description Methods for Time-Weighted DTW objects 
#' 
#' @param x an R object. See \code{\link[dtwSat]{twdtw}}
#' @docType methods
#' @examples
#' alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)
#' is.twdtw(alig)
#' @export
is.twdtw = function(x){
  return(inherits(x,"twdtw"))
}

#' @title Method for Time-Weighted DTW 
#' 
#' @description Methods for Time-Weighted DTW objects 
#' 
#' @param x a \code{\link[dtwSat]{twdtw}} object 
#' @param ... additional arguments passed to print 
#' 
#' @docType methods
#' @examples
#' alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)
#' print(alig)
#' @export
print.twdtw = function(x,...){
  head = "Time-Weighted DTW alignment object\n"
  size = sprintf("Alignment size (query x template): %d x %d\n", x$N, x$M)
  alignments = sprintf("Number of alignment: %d\n", nrow(x$alignments))
  call = sprintf("Call: %s\n", deparse(x$call))
  cat(head,size,alignments,call)
}


#' @title DTW backtrack
#' 
#' @description This function preforms the backtrack starting from  
#' a given index of the last line in the global cost matrix. 
#' 
#' @param alignment A twdtw or dtw alignment object. 
#' @param jmin An integer. The index of the last line in the 
#' global cost matrix. 
#' @docType methods
#' @export
.kthbacktrack = function(alignment, jmin=NULL) {
  
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
  
  out = list()
  out$index1 = ii
  out$index2 = jj
  return(out)
}







