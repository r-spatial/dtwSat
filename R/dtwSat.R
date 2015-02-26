
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



#' @title Satellite image time series analysis. 
#' 
#' @description This function applies an open boundary DTW analysis and 
#' retrieves all possible alignments of a query within a template.
#' It also allows to compute some statistics for each dtw match.
#' 
#' @param query A zoo object with the query time series.
#' @param template A zoo object with the template time series. 
#' It must be larger than the query.
#' @param theta A real between 1 and 0. It is a parameter for dtw local
#' cost matrix computation. For theta equal 0 the function computes a 
#' normal dtw, for greater values it includes the time cost in the local 
#' matrix. Default is 0.
#' @param span A real between 1 and 0. The span says how much the length 
#' of a dtw match can be diferent from the length of the pattern. Default 
#' is 2/3.
#' @param normalize A logical, default is TRUE. Setting TRUE it divides 
#' the DTW distance by the pattern length.
#' @param satStat A logical, default is FALSE. Setting TRUE the function 
#' also computes some statistics for each dtw match.
#' @docType methods
#' @export
timeSeriesAnalysis = function(query, template, theta=0, span=2/3,
                              normalize=TRUE, satStat=FALSE)
{
  
  if(!is.zoo(query))
    stop("Missing zoo object. The query must be a zoo object.")
  if(!is.zoo(template))
    stop("Missing zoo object. The template must be a zoo object.")
  
  if( span < 0 & 1 < span )
    stop("Error: span must be a number between 0 and 1.")
  
  if( theta < 0 & 1 < theta )
    stop("Error: theta must be a number between 0 and 1.")
  
  tx = index(query)
  x  = as.numeric(query)
  ty = index(template)
  y  = as.numeric(template)
  patternLength = abs(tx[length(tx)] - tx[1])
  
  # Step 1. Compute the open boundary DTW between the query and the template
  alignment = .dtwSat(query, template, theta, step.matrix = symmetric1, window.function = noWindow) 

  # Step 2. Retrieve the end point of each path (min points in the last line of the cost matrix) 
  d = alignment$costMatrix[alignment$N,1:alignment$M]
  NonNA = which(!is.na(d))
  diffd = diff(d[NonNA])
  minPoints = NonNA[which(diffd[-length(diffd)] < 0 & diffd[-1] >= 0)] + 1
  
#   gp = ggplot(data = data.frame(x=ty[NonNA], y=d[NonNA]), aes( x = as.Date(x), y = y )) + 
#     ylim(c(0,10)) + 
#     geom_line(size = .5, colour="black", fill="black") + 
#     geom_point(data=data.frame(x=as.Date(ty[b]), y=d[b]), color="red") + 
#     scale_x_date(labels = date_format("%Y"), breaks = date_breaks("year")) +
#     xlab("Time") + 
#     ylab("DTW distance") +
#     ggtitle(paste("Pattern -",patternName)) + 
#     theme(text=element_text(size = 10, family="Helvetica"), plot.title = element_text(size = 10),
#           axis.text.x = element_text(size = 10, family="Helvetica"), axis.text.y = element_text(size = 10, family="Helvetica"))
#   print(gp)
#   
  
  out = do.call("rbind",lapply(minPoints, function(b){
    # Step 3. Map whole possible min paths (k-th paths)
    map = kthbacktrack(alignment, b)
    # Step 4. Retriave the starts of each path
    a = map$index2[1]
    # Step 5. Remove tiny matches 
    lengthDays = abs(ty[b] - ty[a])
    validSubsec = (1 - span) * patternLength <= lengthDays & lengthDays <= (1+span) * patternLength
    if(!validSubsec)
      return(data.frame(dtw.a = NA, dtw.b = NA, dtw.from = as.Date(NA), dtw.to = as.Date(NA),
                        dtw.cost = NA, dtw.timeCost = NA, stringsAsFactors = NA))
    # Step 6. Compute the time cost of each path
    t1 = as.numeric(format(tx[map$index1], "%j"))
    t2 = as.numeric(format(ty[map$index2], "%j"))
    timeCost = theta * sum(abs(t1 - t2)) / 366
    dtwDist = d[b] - timeCost
    if(normalize)
      dtwDist = dtwDist / length(tx)
    
    return(data.frame(dtw.a = a, dtw.b = b, dtw.from = as.Date(ty[a]), dtw.to = as.Date(ty[b]),
                      dtw.cost = dtwDist, dtw.timeCost = timeCost, stringsAsFactors = FALSE))
    
  })) # End mapping loop

  # Remove null lines
  nonnull = !apply(out, 1, is.na)[1,]
  return(out[nonnull,])
}

#' @title Satellite image time series analysis. 
#' 
#' @description This function applies an open boundary DTW analysis and 
#' retrieves all possible alignments of a query within a template.
#' It also allows to compute some statistics for each dtw match.
#' 
#' @param query A zoo object with the query time series.
#' @param template A zoo object with the template time series. 
#' It must be larger than the query.
#' @param theta A real between 1 and 0. It is a parameter for dtw local
#' cost matrix computation. For theta equal 0 the function computes a 
#' normal dtw, for greater values it includes the time cost in the local 
#' matrix. Default is 0.
#' @param span A real between 1 and 0. The span says how much the length 
#' of a dtw match can be diferent from the length of the pattern. Default 
#' is 2/3.
#' @param normalize A logical, default is TRUE. Setting TRUE it divides 
#' the DTW distance by the pattern length.
#' @docType methods
#' @export
timeSeriesAnalysis2 = function(query, template, theta=0, span=2/3,
                              normalize=TRUE)
{
  
  if(!is.zoo(query))
    stop("Missing zoo object. The query must be a zoo object.")
  if(!is.zoo(template))
    stop("Missing zoo object. The template must be a zoo object.")
  
  if( span < 0 & 1 < span )
    stop("Error: span must be a number between 0 and 1.")
  
  if( theta < 0 & 1 < theta )
    stop("Error: theta must be a number between 0 and 1.")
  
  tx = index(query)
  x  = as.numeric(query)
  ty = index(template)
  y  = as.numeric(template)
  
  # Step 1. Compute the open boundary DTW between the query and the template
#   lm = .localCostMatrix(query, template, theta)
#   alignment = dtw(x=lm, step.pattern=symmetric0, # New symmetric with normalization N (see dtw package documentation)
#                   keep.internals=TRUE,open.begin=TRUE,open.end=TRUE)
  
  alignment = .dtwSat(query, template, theta, step.matrix = symmetric1, window.function = noWindow) 

  
  # Step 2. Retrieve the end point of each path (min points in the last line of the cost matrix) 
  d = alignment$costMatrix[alignment$N,1:alignment$M]
  NonNA = which(!is.na(d))
  diffd = diff(d[NonNA])
  b = NonNA[which(diffd[-length(diffd)] < 0 & diffd[-1] >= 0)] + 1
  
  # Step 3. Map whole possible min paths (k-th paths)
  mapping = lapply(b, function(k){
    alignment$jmin = k
    return(kthbacktrack(alignment))
  }) # End mapping loop
  
  # Step 4. Retriave the starts of each path
  a = unlist(lapply(mapping, function(map){
    return(map$index2[1])
  })) # End a loop
  
  # Step 5. Compute the time cost of each path
  timeCost = unlist(lapply(mapping, function(map){
    t1 = as.numeric(format(tx[map$index1], "%j"))
    t2 = as.numeric(format(ty[map$index2], "%j"))
    return(theta * sum(abs(t1 - t2)) / 366)
  })) # End a loop
  
  # Step 6. Remove tiny matches 
  patternLength = abs(tx[length(tx)] - tx[1])
  lengthDays = abs(ty[b] - ty[a])
  validSubsec = (1 - span) * patternLength <= lengthDays & lengthDays <= (1+span) * patternLength
  a = a[validSubsec]
  b = b[validSubsec]
  timeCost = timeCost[validSubsec]
  dtwDist = d[b] - timeCost
  if(normalize)
    dtwDist = dtwDist / length(tx)
  
    return(data.frame(dtw.a = a, dtw.b = b, dtw.from = as.Date(ty[a]), dtw.to = as.Date(ty[b]),
                      dtw.cost = dtwDist, dtw.timeCost = timeCost, stringsAsFactors = FALSE))
    
}

.statTimeSat = function(y, ty, a, b){

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


.dtwSat =  function (query, template, theta, step.matrix = symmetric1, window.function = noWindow) 
{
  tx = index(query)
  x  = as.numeric(query)
  ty = index(template)
  y  = as.numeric(template)
  tx = as.numeric(format(tx, "%j")) / 366
  ty = as.numeric(format(ty, "%j")) / 366
  lm = proxy::dist((1-theta)*x,(1-theta)*y,method="euclidean")
  tm = proxy::dist(theta*tx,theta*ty,method="euclidean")
  lm = lm + tm
  
  lm <- rbind(0, lm)
  n = nrow(lm)
  m = ncol(lm)
  cm = matrix(NA, nrow=n, ncol=m)
  cm[1,] = 0    
  wm = matrix(FALSE, nrow = n, ncol = m)
  wm[window.function(row(wm), col(wm), query.size = n, reference.size = m)] = TRUE
  out = .Call("computeCM_Call", PACKAGE="dtw", wm, lm, cm, step.matrix)
  out$costMatrix = out$costMatrix[-1,]
  out$directionMatrix = out$directionMatrix[-1,]
  out$stepPattern = step.matrix
  out$N = n-1
  out$M = m
  class(out) <- "dtw"
  return(out)
}


.localCostMatrix = function(query, template, theta){
  tx = index(query)
  x  = as.numeric(query)
  ty = index(template)
  y  = as.numeric(template)
  tx = as.numeric(format(tx, "%j")) / 366
  ty = as.numeric(format(ty, "%j")) / 366
  lm = proxy::dist((1-theta)*x,(1-theta)*y,method="euclidean")
  tm = proxy::dist(theta*tx,theta*ty,method="euclidean")
  lm = lm + tm
  return(lm)
}

.globalCostMatrix =  function (lm, step.matrix = symmetric1, window.function = noWindow) 
{
  lm <- rbind(0, lm)
  n = nrow(lm)
  m = ncol(lm)
  cm = matrix(NA, nrow=n, ncol=m)
  cm[1,] = 0    
  wm = matrix(FALSE, nrow = n, ncol = m)
  wm[window.function(row(wm), col(wm), query.size = n, reference.size = m)] = TRUE
  out = .Call("computeCM_Call", PACKAGE="dtw", wm, lm, cm, step.matrix)
  out$costMatrix = out$costMatrix[-1,]
  out$directionMatrix = out$directionMatrix[-1,]
  out$stepPattern = step.matrix
  out$N = n-1
  out$M = m
  class(out) <- "dtw"
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
        subsequence = aux[k,grep(names(aux), pattern="dtw.cost")]
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
          dtwcost = subsequence[grep(names(subsequence), pattern="dtw.cost")]
          dtwcost[is.na(dtwcost)] = THRESHOLD 
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
          I = which(as.numeric(format(x$dtw.to, "%Y"))==y)
          if(length(I)==0)
            return(NA)
          return(min(x$dtw.cost[I]))
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
  names(x) = paste("dtw.cost", seq_along(x), sep=".")
  
  return(data.frame(classNames, x))
}


.lowestDTWCostInOverlapping = function(x, year, startDate, endDate, overlapping)
{
    I = which( year == as.numeric(format(x$dtw.to,"%Y")) )
    
    # Test minimum overlapping 
    flag = unlist(lapply(I, function(i){
      a = x$dtw.from[i]
      if(a < startDate)
        a = startDate
      b = x$dtw.to[i]
      if(endDate < b)
        b = endDate
      
      if(a >= endDate | b <= startDate)
        return(FALSE)
      
      subsequencesOverlapping = as.numeric(b - a) / as.numeric(endDate - startDate)
      return(subsequencesOverlapping >= overlapping)
    }))
    
    if(!any(flag))
      return(data.frame(year = year, dtw.cost = NA, dtw.timeCost = NA, 
                        stringsAsFactors = FALSE))
    
    # Return the lowest DTW cost
    I = I[which.min(x$dtw.cost[I[flag]])]
    return(data.frame(year = year, dtw.cost = x$dtw.cost[I], 
                      dtw.timeCost = x$dtw.timeCost[I], stringsAsFactors = FALSE))

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








#' @title Satellite image time series analysis. 
#' 
#' @description This function applies an open boundary DTW analysis and 
#' retrieves all possible alignments of a query within a template.
#' It also allows to compute some statistics for each dtw match.
#' 
#' @param query A zoo object with the query time series.
#' @param template A zoo object with the template time series. 
#' It must be larger than the query.
#' @param theta A real between 1 and 0. It is a parameter for dtw local
#' cost matrix computation. For theta equal 0 the function computes a 
#' normal dtw, for greater values it includes the time cost in the local 
#' matrix. Default is 0.
#' @param span A real between 1 and 0. The span says how much the length 
#' of a dtw match can be diferent from the length of the pattern. Default 
#' is 2/3.
#' @param normalize A logical, default is TRUE. Setting TRUE it divides 
#' the DTW distance by the pattern length.
#' @param satStat A logical, default is FALSE. Setting TRUE the function 
#' also computes some statistics for each dtw match.
#' @docType methods
#' @export
timeSeriesAnalysis2 = function(query, template, theta=0, span=2/3,
                              normalize=TRUE, satStat=FALSE)
{
  
  if(!is.zoo(query))
    stop("Missing zoo object. The query must be a zoo object.")
  if(!is.zoo(template))
    stop("Missing zoo object. The template must be a zoo object.")
  
  if( span < 0 & 1 < span )
    stop("Error: span must be a number between 0 and 1.")
  
  if( theta < 0 & 1 < theta )
    stop("Error: theta must be a number between 0 and 1.")
  
  tx = index(query)
  x  = as.numeric(query)
  ty = index(template)
  y  = as.numeric(template)
  
  # Step 1. Compute the open boundary DTW between the query and the template
  lm = .localCostMatrix(query, template, theta)
  alignment = dtw(x=lm, step.pattern=symmetric0, # New symmetric with normalization N (see dtw package documentation)
                  keep.internals=TRUE,open.begin=TRUE,open.end=TRUE)
  
  # Step 2. Retrieve the end point of each path (min points in the last line of the cost matrix) 
  d = alignment$costMatrix[alignment$N,1:alignment$M]
  NonNA = which(!is.na(d))
  diffd = diff(d[NonNA])
  b = NonNA[which(diffd[-length(diffd)] < 0 & diffd[-1] >= 0)] + 1
  
  #   gp = ggplot(data = data.frame(x=ty[NonNA], y=d[NonNA]), aes( x = as.Date(x), y = y )) + 
  #     ylim(c(0,10)) + 
  #     geom_line(size = .5, colour="black", fill="black") + 
  #     geom_point(data=data.frame(x=as.Date(ty[b]), y=d[b]), color="red") + 
  #     scale_x_date(labels = date_format("%Y"), breaks = date_breaks("year")) +
  #     xlab("Time") + 
  #     ylab("DTW distance") +
  #     ggtitle(paste("Pattern -",patternName)) + 
  #     theme(text=element_text(size = 10, family="Helvetica"), plot.title = element_text(size = 10),
  #           axis.text.x = element_text(size = 10, family="Helvetica"), axis.text.y = element_text(size = 10, family="Helvetica"))
  #   print(gp)
  #   
  
  # Step 3. Map whole possible min paths (k-th paths)
  mapping = lapply(b, function(k){
    alignment$jmin = k
    return(kthbacktrack(alignment))
  }) # End mapping loop
  
  # Step 4. Retriave the starts of each path
  a = unlist(lapply(mapping, function(map){
    return(map$index2[1])
  })) # End a loop
  
  # Step 5. Compute the time cost of each path
  timeCost = unlist(lapply(mapping, function(map){
    t1 = as.numeric(format(tx[map$index1], "%j"))
    t2 = as.numeric(format(ty[map$index2], "%j"))
    return(theta * sum(abs(t1 - t2)) / 365)
  })) # End a loop
  
  # Step 6. Remove tiny matches 
  patternLength = abs(tx[length(tx)] - tx[1])
  lengthDays = abs(ty[b] - ty[a])
  validSubsec = (1 - span) * patternLength <= lengthDays & lengthDays <= (1+span) * patternLength
  a = a[validSubsec]
  b = b[validSubsec]
  timeCost = timeCost[validSubsec]
  dtwDist = d[b] - timeCost
  if(normalize)
    dtwDist = dtwDist / length(tx)
  
  if(!satStat){
    return(data.frame(dtw.a = a, dtw.b = b, dtw.from = as.Date(ty[a]), dtw.to = as.Date(ty[b]),
                      dtw.cost = dtwDist, dtw.timeCost = timeCost, stringsAsFactors = FALSE))
  }
  
  return(data.frame(dtw.a = a, dtw.b = b, dtw.from = as.Date(ty[a]), dtw.to = as.Date(ty[b]),
                    dtw.cost = dtwDist, dtw.timeCost = timeCost, .statTimeSat(y, ty, a, b), stringsAsFactors = FALSE))
  
  
}


