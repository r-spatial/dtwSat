
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
#' cost matrix computation. For theta equal 1 the function computes a 
#' normal dtw, for lower values it includes the time cost in the local 
#' matrix. Default is 1.
#' @param span A real between 1 and 0. The span says how much the length 
#' of a dtw match can be diferent from the length of the pattern. Default 
#' is 2/3.
#' @param normalize A logical, default is TRUE. Setting TRUE it divides 
#' the DTW distance by the pattern length.
#' @param satStat A logical, default is FALSE. Setting TRUE the function 
#' also computes some statistics for each dtw match.
#' @docType methods
#' @export
timeSeriesAnalysis = function(query, template, theta=1.0, span=2/3,
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
  alignment = dtw(x=x, y=y, tx=tx, ty=ty, theta=theta,
                  step=symmetric0, # New symmetric with normalization N (see dtw package documentation)
                  keep=TRUE,open.begin=TRUE,open.end=TRUE)
  
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
    return((1 - theta) * sum(abs(t1 - t2)) / 365)
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






#' @title Satellite image time series classification. 
#' 
#' @description This function classify the time segments based on the 
#' open boundary DTW analysis results (i.e. results of timeSeriesAnalysis 
#' function).
#' 
#' @param dtwResults A data frame with the based on the open boundary DTW
#' analysis results. It must be larger than the query.
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
#' @docType methods
#' @export
timeSeriesClassifier = function(dtwResults, from, to, by=12, 
                                overlapping=0.4, threshold=4.0)
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
    lowesCous = lapply( dtwResults, function(x){
      out = do.call("rbind", lapply(seq_along(years), function(k){
        .lowestDTWCostInOverlapping(x, years[k], startDates[k], endDates[k], overlapping)
      }))
    })
    
#     # Update probabilit
#     lowesCous = lapply( years, function(y){
#       lapply( dtwResults, function(x) x$dtw.cost[x$year==y })
#     })
#     
#     dtwYX = distMat[1,2]
#     dtwX = lowesCous[[1]]$dtw.cost[1]
#     dtwY = lowesCous[[1]]$dtw.cost[2]
#     dtwXY = (dtwYX * dtwX) / dtwY
#     

    # Find the best classification based on DTW cost
    lowesCous = data.frame(lowesCous)
    bestClass = do.call("rbind", lapply(seq_along(years), function(k){
        subsequence = lowesCous[k,grep(names(lowesCous), pattern=".dtw.cost")]
        out = .bestDTWClass(subsequence, threshold=threshold)
        data.frame(year=years[k], out)
    }))

    return(bestClass)
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






# TODO: Add comment
# 
# Author: Maus, Victor
###############################################################################
#
#
#
#     ############### EXTRACTED FROM DTW PACKAGE - By Victor Maus ##################
#            These functions are needed to compute the k-th minimal paths
#
#
#
#
#
## Extract rows belonging to pattern no. sn
## with first element stripped
## in reverse order

.extractpattern <- function(sp,sn) {
  sbs<-sp[,1]==sn;	# pick only rows beginning by sn
  spl<-sp[sbs,-1,drop=FALSE];
  # of those: take only column Di, Dj, cost
  # (drop first - pattern no. column)
  
  nr<-nrow(spl);	# how many are left
  spl<-spl[nr:1,,drop=FALSE];	# invert row order
  
  return(spl);
}


###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino@gmail.com>               #
#       Laboratory for Biomedical Informatics                 #
#       University of Pavia - Italy                           #
#       www.labmedinfo.org                                    #
#                                                             #
#   $Id: backtrack.R 168 2008-07-11 05:52:05Z tonig $		  #
#                                                             #
###############################################################


########################################
## Backtrack the steps taken - internal

`kthbacktrack` <- function(gcm) {
  
  dir<-gcm$stepPattern;
  npat <- attr(dir,"npat");
  
  n <- nrow(gcm$costMatrix);
  m <- ncol(gcm$costMatrix);
  
  i <- n;
  j <- gcm$jmin
  
  
  
  ## drop rows with (0,0) deltas 
  nullrows <- dir[,2]==0 & dir[,3] ==0 ;
  tmp <- dir[!nullrows,];
  
  ## Pre-compute steps
  stepsCache <- list();  
  for(k in 1:npat) {
    stepsCache[[k]] <- .extractpattern(tmp,k);
  }
  
  
  ## mapping lists
  ii<-c(i);
  jj<-c(j);
  
  
  repeat {
    ## cross fingers for termination
    if(i==1){# && j==1) {
      break; 	
    }
    ## direction taken
    s<-gcm$directionMatrix[i,j];
    if(is.na(s)) {
      break;
    }
    
    ## undo the steps
    
    steps<-stepsCache[[s]];
    ns<-nrow(steps);
    
    ## In some rare cases (eg symmetricP0), ns will be 1
    ## R indexing rules make k==0 a no-op anyway
    for(k in 1:ns) {
      ## take note of current cell, prepending to mapping lists
      if(i-steps[k,1] > 0) {       # Modified from original function
        ii <- c(i-steps[k,1],ii);  # Modified from original function
        jj <- c(j-steps[k,2],jj);  # Modified from original function
      }                            # Modified from original function
      ## All sub-steps are visited & appended; we have dropped (0,0) deltas
    }
    
    ## And don't forget where we arrived to
    i <- ii[1]#i-steps[ns,1]; # Modified from original function
    j <- jj[1]#j-steps[ns,2]; # Modified from original function
  }
  ############################################################################
  
  out<-list();
  out$index1<-ii;
  out$index2<-jj;
  
  return(out);
}
