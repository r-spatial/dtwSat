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
#' @param timeline A time line vector for the internalsput.
#' @param frequency A numeric with the frequency for internalsput time series.
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
#' @examples
#' ##
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




#' @title Classification accuracy  
#' 
#' @description This function retrieves the accuracy of classification
#' based on the a given dataset of controll points.
#' 
#' @param predicted A vector of prediction. 
#' @param reference A vector of reference for the predicted classes. 
#' @docType methods
#' @examples
#' ##
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
