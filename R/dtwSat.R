
#' @title Satellite image time series smoothing. 
#' 
#' @description This function applies a smoothing algorithm to 
#' the satellite image time series. It allows to compute the 
#' discreat wavelet smoothing for different levels and the 
#' Savitzky–Golay smoothing. This function only works for constant 
#' frequency data, then informa a constant frequency time series.
#' If it is not the case, either informe a constant frequency time line 
#' or a frequency. 
#' 
#' @param x Either a time line vector or a zoo object with the 
#' time series.
#' @param y A numeric vector.
#' @param timeline A time line vector for the output.
#' @param frequency A numeric with the frequency for output time series.
#' @param method A character vector for the smoothing methods, either 
#' wavelet and level, eg. c("wavelet","1") or Savitzky–Golay, polynomial 
#' and a window size, e.g. c("sg", "2", "3").
#' @docType methods
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
  }else if(method[1]=="sg") {
    stop("I'm sorry. It is not implemented yet!")
  } else {
    stop("Missing smoothing method. Please choose between discrete wavelet and Savitzky–Gola!")
  }
  
  return(zoo(sy, sty))

}


#' @title Multidimensional Time-Weighted DTW analysis. 
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves one or more possible alignments of a query within 
#' a time series.
#' 
#' @param query A zoo object with the multidimensional time series.
#' @param template A zoo object with the template time series. 
#' It must be iguel or be equal or longer than the length of the query and 
#' the same number of dimensions.
#' @param weight A character. "linear" for linear weight or "logistic" 
#' for logistic weight. Default is NULL that runs the original dtw method.
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is "Euclidean" See \code{\link[proxy]{dist}} in package 
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
    endPoints = sort(endPoints)
    if(is.null(alignments))
      alignments = length(endPoints)
    if(length(endPoints) > alignments)
      endPoints = endPoints[1:alignments]
    # Map low cost paths (k-th paths)
    mapping = lapply(endPoints, function(b){
      return(kthbacktrack(out, b))
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
  out$call = match.call()
#   out$stepPattern = 
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
#' @param x twdtw object
#' @param n An integer. The alignment. Default is 1
#' @param d An integer. The dimension for point-by-point 
#' comparison, with matching lines \code{\link{dtwPlotThreeWay}} 
#' in package \pkg{dtw}. Default is 1
#' @param type see \code{\link{plot.dtw}} in package \pkg{dtw}. 
#' Default is "density"
#' @param ... see \code{plot.dtw} in package \pkg{dtw}
#' @docType methods
#' @export
plotTwdtw = function(x, n=1, d=1, type="density", ...){
  if(0 > n & n > length(x$mapping))
    stop(paste("n must be an integer between 1 and the number of alignments:",length(x$mapping)))
  if(0 > d & d > ncol(x$query) )
    stop(paste("d must be an integer between 1 and the number of cols in the query:",ncol(x$query)))
  obj = .dtwObj(x, n=1, d=1)
  dtw::dtwPlot(obj, type="density", ...)
}

# PLOT ALL ALIGNMENTS IN THE SAME TERREAIN
# RUN LIST OF PATTERNS AND RETRIENVE A DATA.FRAME

.dtwObj = function(x, n=1, d=1){
  if(class(x)!="twdtw")
    stop("Missing a twdtw object.")
  res = x
  res$jmin      = x$alignments$to
  res$index1    = x$mapping[[n]]$index1
  res$index2    = x$mapping[[n]]$index2
  res$query     = x$query[,d]
  res$reference = x$template[,d]
  return(res)
}

#' @title Function to compute phenological metrics. 
#' 
#' @description This compute phenological metrics from a time series.
#' 
#' @param y A zoo object with the query time series.
#' @param a A vector. Indices in y of alignment starts.
#' @param b A vector. Indices in y of alignment ends.
#' @docType methods
#' @export
statTimeSat = function(y, a=NULL, b=NULL){
  if(!is.zoo(y))
    stop("Missing zoo object. The query must be a zoo object.")
  ty = as.Date(index(y))
  y  = as.numeric(y)
  if(is.null(a)) a = 1
  if(is.null(b)) b = length(y)
  out = do.call("rbind", lapply(seq_along(a), function(k){
    imax = which.max(y[a[k]:b[k]]) + a[k] - 1
    dy = diff(y[a[k]:b[k]])
    ipics = which(dy[-1] < 0 & 0 < dy[-length(dy)]) + a[k] + 1
    if( length(ipics)==0 )
      ipics = imax
    iminL = which.min(y[a[k]:ipics[1]]) + a[k] - 1
    iminR = which.min(y[ipics[length(ipics)]:b[k]]) + ipics[length(ipics)] - 1
    dateLeftMin = ty[iminL]
    dateRightMin = ty[iminR]
    dateMax = ty[imax]
    pics = length(ipics)
    slength = as.numeric(ty[b[k]] - ty[a[k]])
    leftMin = y[iminL]
    rightMin = y[iminR]
    smax = y[imax]
    smin = mean(leftMin, rightMin, na.rm=TRUE)
    smplitude = smax - smin
    smean = mean(y[a[k]:b[k]], na.rm=TRUE)
    ssd = sd(y[a[k]:b[k]], na.rm=TRUE)
    leftSlop = (y[ipics[1]] - leftMin) / abs(as.numeric(ty[ipics[1]] - ty[iminL]))
    rightSlop = (rightMin - y[ipics[length(ipics)]]) / abs(as.numeric(ty[ipics[length(ipics)]] - ty[iminL]))
    totInt = sum(y[a[k]:b[k]], na.rm=TRUE)
    seasonInt = sum(y[a[k]:b[k]]-smin, na.rm=TRUE)
    return(
      data.frame(
        dateLeftMin = dateLeftMin,
        dateRightMin = dateRightMin,
        dateMax = dateMax,
        pics = pics,
        length = slength,
        leftMin = leftMin,
        rightMin = rightMin,
        max = smax,
        min = smin,
        mean = smean,
        sd = ssd,
        leftSlop = leftSlop,
        rightSlop = rightSlop,
        totInt = totInt,
        seasonInt = seasonInt
      ))
  }))
  
  return(out)
}


#' @title Satellite image time series classification. 
#' 
#' @description This function classify the time segments based on the 
#' open boundary DTW analysis results (i.e. results of timeSeriesAnalysis 
#' function).
#' 
#' @param dtwResults A list of dtw . 
#' @param from A character with the start date for the finalclassification.
#' The format is yyyy-mm-dd.
#' @param to A character with the end date for the finalclassification. The 
#' format is yyyy-mm-dd.
#' @param by An integer. The length in months of the time subsequences
#' for the final classification. Default is 12 months.
#' @param overlapping A real between 0 and 1 representing the percentage of
#' minimum overlapping between the DTW mach and the subsequences for the 
#' final classification. Default is 40\%.
#' @param threshold A real with the DTW threshold, i.e. the maximum DTW cost 
#' for consideration.
#' @param sortBydtw A logical. If TRUE the function returns the best 
#' classification for each period of time. If FALSE the functions 
#' returns the lowest dtw for each class. Default is TRUE.
#' @param aggregateByClass A logical. Aggregate results by common class. 
#' Default is TRUE
#' @docType methods
#' @export
timeSeriesClassifier = function(dtwResults, from, to, by=12,
                                overlapping=0.4, threshold=4.0, 
                                sortBydtw=TRUE, aggregateByClass=TRUE)
{
    if(!is.list(dtwResults))
      stop("Missing a list. The parameter dtwResults must be a list.")
    
    if( !(is.character(from) | class(from)=="Date") )
      stop("Missing a character. The parameter from must be a character or a Date in the formt yyyy-mm-dd.")
    
    if(!(is.character(to) | class(to)=="Date") )
      stop("Missing a character. The parameter to must be a character or a Date in the formt yyyy-mm-dd.")
    
    if( overlapping < 0 & 1 < overlapping )
      stop("Error: overlapping must be a number between 0 and 1.")
    
    from = as.Date(from, origin="1970-01-01")
    to = as.Date(to, origin="1970-01-01")
    startDates = seq(from, to, by = paste(by, "month"))
    endDates = seq(startDates[2], to+by-1, by = paste(by, "month"))
    years = as.numeric(format(endDates,"%Y"))
    
    className = names(dtwResults)
    if(is.null(className))
      className = paste("class",seq_along(dtwResults),sep="_")
    
    # Find the lowest cost for each class
    aux = lapply( dtwResults, function(x){
       out = do.call("rbind", lapply(seq_along(years), function(k){
         .lowestDTWCostInOverlapping(x, years[k], startDates[k], endDates[k], overlapping)
       }))
    })
    aux = data.frame(aux)
    
    if(sortBydtw){
      # Find the best classification based on DTW cost
      res = do.call("rbind", lapply(seq_along(years), function(k){
        subsequence = aux[k,grep(names(aux), pattern="distance")]
        out = .bestDTWClass(subsequence, threshold=threshold)
        data.frame(year=years[k], out)
      }))
      return(res)
    }

    
    if(aggregateByClass){
      res = lapply(c("years", unique(className)), function(name){
        if(name=="years")
          return(years)
        out = unlist(lapply(seq_along(years), function(k){
          subsequence = aux[k,grep(names(aux), pattern=name)]
          dtwcost = subsequence[grep(names(subsequence), pattern="distance")]
          dtwcost[is.na(dtwcost)] = threshold 
          return(min(as.numeric(dtwcost), na.rm=TRUE))
        }))
        return(out)
      })
      names(res) = c("year", unique(className))
      return(res)
    }
    
    names(dtwResults) = paste(names(dtwResults),seq_along(dtwResults),sep="_")
    
    res = lapply(dtwResults, function(x){
      unlist(lapply(years, function(y){
          I = which(as.numeric(format(x$to, "%Y"))==y)
          if(length(I)==0)
            return(NA)
          return(min(x$distance[I]))
      }))
    })
    
    return(data.frame(year=years, res))  

}


.bestDTWClass = function(x, threshold)
{
  x[x > threshold] = NA
  x = x[order(x)]
  
  classNames = data.frame(lapply(names(x), function(name){
    unlist(strsplit(name, split="\\."))[1]
  }), stringsAsFactors = FALSE )
  
  classNames[is.na(x)] = "other"
  
  names(classNames) = paste("class", seq_along(x), sep=".")
  names(x) = paste("distance", seq_along(x), sep=".")
  
  return(data.frame(classNames, x))
}


.lowestDTWCostInOverlapping = function(x, year, startDate, endDate, overlapping)
{
    I = which( year == as.numeric(format(x$to,"%Y")) )
    
    # Test minimum overlapping 
    flag = unlist(lapply(I, function(i){
      a = x$from[i]
      if(a < startDate)
        a = startDate
      b = x$to[i]
      if(endDate < b)
        b = endDate
      
      if(a >= endDate | b <= startDate)
        return(FALSE)
      
      subsequencesOverlapping = as.numeric(b - a) / as.numeric(endDate - startDate)
      return(subsequencesOverlapping >= overlapping)
    }))
    
    if(!any(flag))
      return(data.frame(year = year, 
                        distance = NA,
                        stringsAsFactors = FALSE))
    
    # Return the lowest DTW cost
    I = I[which.min(x$distance[I[flag]])]
    return(data.frame(year = year, 
                      distance = x$distance[I], 
                      stringsAsFactors = FALSE))

}


#' @title Computes the accuracy of a classification
#' 
#' @description This function computs the accuracy a classification based
#' on the a given dataset of controll points.
#' 
#' @param predicted A vector of prediction. 
#' @param reference A vector of reference for the predicted classes. 
#' @docType methods
#' @export
computeAccuracy = function(predicted, reference)
{

  if( length(predicted)!=length(reference) )
    stop("The vectors 'predicted and 'reference' must have the same length.")
  
  factorLevels = unique( c(reference,predicted) )
  
  # Global accuracy
  global_accuracy = sum(predicted==reference) / length(reference)
  
  a = factor(predicted, levels = factorLevels)
  b = factor(reference, levels = factorLevels)
    
  # Check predicted classes against reference class labels
  confusion_table = table(predicted=a, reference=b)
  
  return(list(global_accuracy = global_accuracy, 
              confusion_table = confusion_table))
}



#' @title Computes the backtrack of dtw
#' 
#' @description This function computs the backtrack starting from  
#' a given index of the last line of the global cost matrix. 
#' 
#' @param alignment A dtw alignment object. 
#' @param jmin An integer. index of the last line of the global 
#' cost matrix. 
#' @docType methods
#' @export

kthbacktrack = function(alignment, jmin=NULL) {
  
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







