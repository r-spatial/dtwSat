#' @title Wavelet filter
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a smoothing algorithm to
#' the time series. It computes the a discreat wavelet 
#' smoothing for each dimension in the imput time series.
#' 
#' @param x A zoo object with the time series
#' @param timeline A time line vector for the internalsput
#' @param frequency A numeric with the frequency for internalsput time series
#' @param wf Name of the wavelet filter to use in the decomposition. 
#' Default is "la8"
#' @param J Specifies the depth of the decomposition. This must be a number 
#' less than or equal to log(length(x),2). Default is 1
#' @param boundary Character string specifying the boundary condition. 
#' Default is "periodic"
#' @param ... see parameters of \code{\link[waveslim]{mra}} in the 
#' packege \pkg{waveslim}
#' @docType methods
#' @return object of class \code{\link[zoo]{zoo}} 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' ## Plot raw EVI and filtered EVI
#' # require(ggplot2)
#' #df = data.frame(Time=index(template), value=template$evi, variable="Raw")
#' #df = rbind( df, data.frame(Time=index(sy), value=sy$evi, variable="Wavelet filter") )
#' #ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
#' #        geom_line() +
#' #       ylab("EVI")
#'     
#' ## Plot all filter bands
#' # require(ggplot2)
#' # require(reshape2)
#' #df = melt(data.frame(Time=index(sy), sy), id="Time")
#' #ggplot(df, aes(x=Time, y=value, group=variable, colour=variable)) +
#' #       geom_line() 
#'   
#' @export
waveletSmoothing = function(x, timeline=NULL, frequency=NULL, 
                            wf = "la8", J=1, boundary = "periodic", ...)
{
  if( missing(x) )
    stop("Missing either a numeric vector or a zoo object.")
  if(!is(x, "zoo"))
    stop("x is not a zoo object")
  
  if(is.null(timeline)){
    if(is.null(frequency)){
      timeline = index(x)
    }else{
      timeline = seq(min(index(x)), max(index(x)), by=frequency)
    }
  }
  
  # Linear interpolation of gaps
  I = which(!is.na(timeline) & !duplicated(timeline))
  timeline = timeline[I]
  xx = zoo(, timeline)
  xx = merge(x, xx)
  xx = na.approx(xx)
  xx = xx[timeline,]

  # Smoothing 
  df = lapply(as.list(xx), function(d){
    waveslim::mra(x=as.numeric(d), wf=wf, J=J, boundary=boundary)[[paste0("S",J)]]
  })
  res = zoo(data.frame(df), index(xx))
  res
}


#' @title Satellite image time series classification. 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
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
#' @return object of class \code{\link[base]{data.frame}} 
#' @examples
#' ##
#' 
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


