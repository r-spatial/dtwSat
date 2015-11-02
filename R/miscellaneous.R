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
#### UTILITY FUNCTIONS



#' @title Wavelet filter
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a smoothing algorithm to
#' the time series. It computes a discreat wavelet 
#' smoothing for each dimension in the imput time series.
#' 
#' @param timeseries A \code{\link[zoo]{zoo}} object with the time series
#' @param timeline A vector of dates for the output time series.
#' It must have a regular frequency. 
#' @param frequency The frequncy for the output time series
#' @param wf Name of the wavelet filter used in the decomposition. 
#' Default is "la8"
#' @param J Specifies the depth of the decomposition. This must be a number 
#' less than or equal to log(length(x),2). Default is 1
#' @param boundary Character string specifying the boundary condition. 
#' Default is "periodic". See parameters of \code{\link[waveslim]{mra}}.
#' @param ... see parameters of \code{\link[waveslim]{mra}} in the 
#' packege \pkg{waveslim}
#' @param x is deprecated, please use \code{timeseries} instead
#' 
#' @docType methods
#' @return object of class \code{\link[zoo]{zoo}} 
#' 
#' @seealso \link[waveslim]{mra}
#' 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' plot(template$evi, xlab="Time", ylab="EVI")
#' lines(sy$evi, col="red")
#' 
#' ## Plot raw EVI and filtered EVI
#' #require(ggplot2)
#' #evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
#' #gp = autoplot(evi, facets = NULL) + xlab("Time")
#' #gp
#'     
#' ## Plot all filter bands
#' #require(ggplot2)
#' #gp = autoplot(sy, facets = NULL) + xlab("Time")
#' #gp
#' 
#' @export
waveletSmoothing = function(x=NULL, timeseries=NULL, timeline=NULL, frequency=NULL, 
                            wf = "la8", J=1, boundary = "periodic", ...)
{
  
  if (!missing(x)){
    warning("argument x is deprecated, please use timeseries instead", call. = FALSE)
    timeseries = x
  }
  
  if(!is(timeseries, "zoo"))
    stop("timeseries is not a zoo object")
  
  if(is.null(timeline)){
    if(is.null(frequency)){
      timeline = index(timeseries)
    }else{
      timeline = seq(min(index(timeseries)), max(index(timeseries)), by=frequency)
    }
  }
  
  # Linear interpolation of gaps
  I = which(!is.na(timeline) & !duplicated(timeline))
  timeline = timeline[I]
  y = zoo(order.by = timeline)
  y = merge(timeseries, y)
  y = na.approx(y)
  y = y[timeline,]

  # Smoothing 
  df = lapply(as.list(y), function(d){
    mra(x=as.numeric(d), wf=wf, J=J, boundary=boundary, ...)[[paste0("S",J)]]
  })
  res = zoo(data.frame(df), index(y))
  res
}


#' @title Create time sequence
#' 
#' @description This function creates a sequence of dates for 
#' each year. The sequences start on January 1st of each year.
#' 
#' @param year A vector with the years. Default 
#' is form 2000 to the system time year \code{format(Sys.time(), ''\%Y'')}
#' @param frequency An integer with the frequency in days. Default is 16 days
#' @docType methods
#' 
#' @return vector of class \code{\link[base]{Dates}} 
#' 
#' @seealso \link[dtwSat]{getDatesFromDOY} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' dates = getModisTimeSequence()
#' dates
#' 
#' @export
createTimeSequence = function(year=2000:format(Sys.time(), "%Y"), frequency=16){
  .Deprecated("getModisTimeSequence")
  res = unlist(lapply(year, function(y){
    days = seq(from = as.Date(paste0(y,"-01-01")), to = as.Date(paste0(y,"-12-31")), by = frequency)
  }))
  res = as.Date(res, origin="1970-01-01")
  res
}

#' @title Create time sequence
#' 
#' @description This function creates a sequence of dates for 
#' each year. The sequences start on January 1st of each year.
#' 
#' @param year A vector with the years. Default 
#' is form 2000 to the system time year \code{format(Sys.time(), ''\%Y'')}
#' @param frequency An integer with the frequency in days. Default is 16 days
#' @docType methods
#' 
#' @return vector of class \code{\link[base]{Dates}} 
#' 
#' @seealso \link[dtwSat]{getDatesFromDOY} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' dates = getModisTimeSequence()
#' dates
#' 
#' @export
getModisTimeSequence = function(year=2000:format(Sys.time(), "%Y"), frequency=16){
  res = unlist(lapply(year, function(y){
    days = seq(from = as.Date(paste0(y,"-01-01")), to = as.Date(paste0(y,"-12-31")), by = frequency)
  }))
  res = as.Date(res, origin="1970-01-01")
  res
}

#' @title Get dates from year and day of the year
#' 
#' @description This function retrieves the date corresponding to year 
#' and day of the year
#' 
#' @param year An vector with the years
#' @param doy An vector with the day of the year
#' @docType methods
#' 
#' @return object of class \code{\link[base]{Dates}} 
#' 
#' @seealso \link[dtwSat]{getModisTimeSequence} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' year = c(2000, 2001)
#' doy = c(366, 365)
#' dates = getDatesFromDOY(year, doy)
#' dates
#'
#' @export
getDatesFromDOY = function(year, doy){
  if(length(year)!=length(doy))
    stop("year and doy are not the same length")
  res = as.Date(paste(year, doy), format="%Y %j", origin="1970-01-01")
  res
}


#' @title Get Modis time index from date
#' 
#' @description This function retrieves the nearest time 
#' index to a date
#' 
#' @param date A \code{\link[base]{Dates}} object
#' @param frequency An integer with the frequency in days. Default is 16 days
#' 
#' @docType methods 
#' 
#' @return An integer 
#' 
#' @seealso \link[dtwSat]{getModisTimeSequence} and \link[dtwSat]{getDatesFromDOY}
#' 
#' @examples
#' i = getModisTimeIndex(date=as.Date("2000-01-01"))
#' i
#'
#' @export
getModisTimeIndex = function(date, frequency=16){
  dates = getModisTimeSequence(frequency=frequency)
  res = which.min(abs(dates-date))
  res
}





#' @title Classify time intervals
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the best alignment within each 
#' interval of classification based on the TWDTW distance
#' 
#' @param x A \code{\link[dtwSat]{dtwSat}} object or 
#' a \code{\link[base]{data.frame}} similar to the slot \code{alignments} in 
#' object \code{\link[dtwSat]{dtwSat-class}} 
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param to A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param by A character with the intevals size, \emph{e.g.} ''6 month''
#' @param breaks A vector of class \code{\link[base]{Dates}}
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between the one alignment and the interval of classification. 
#' Default is 1, \emph{i.e.} 100\%
#' @param threshold A number. The TWDTW threshold, i.e. the maximum TWDTW 
#' cost for consideration. Default is \code{Inf}
#' @param normalized Use normalized TWDTW distance. Default is TRUE
#' @docType methods
#' @return object of class \code{\link[base]{data.frame}} with the best alignment 
#' for each interval
#' @examples
#' malig = mtwdtw(query.list, timeseries=template, weight = "logistic", 
#'          normalize=TRUE, query.length=23, alpha = 0.1, beta = 100)
#'          
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' best_class = classifyIntervals(x=malig, from=from, to=to, by = "6 month",
#'              overlap=.3, threshold=Inf)
#' best_class
#' 
#' 
#' @export
classifyIntervals = function(x, from, to, by, breaks=NULL, overlap=.5, threshold=Inf, normalized)
{
  if (!missing(normalized))
    warning("argument normalized is deprecated and is scheduled to be removed in the next version", 
            call. = FALSE)
  
  if(is(x, "dtwSat"))
    x = getAlignments(x)

  if(!is(x, "data.frame"))
    stop("x is not a data.frame or dtwSat class")
  
  if( overlap < 0 & 1 < overlap )
    stop("overlap out of range, it must be a number between 0 and 1")

  if(is.null(breaks))
    breaks = seq(as.Date(from), as.Date(to), by=by)

  res = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    .bestInterval(x, start=breaks[i-1], end=breaks[i], overlap)
  }))
  
  d = res$distance
  I = d>threshold
  if(any(I)){
    res$query[I] = "unclassified"
    res$distance[I] = Inf 
  }
  res
}


#' @title Classify time intervals
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the best alignment within each 
#' interval of classification based on the TWDTW distance
#' 
#' @param x A \code{\link[dtwSat]{dtwSat}} object or 
#' a \code{\link[base]{data.frame}} similar to the slot \code{alignments} in 
#' object \code{\link[dtwSat]{dtwSat-class}} 
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param to A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param by A character with the intevals size, \emph{e.g.} ''6 month''
#' @param breaks A vector of class \code{\link[base]{Dates}}
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between the one alignment and the interval of classification. 
#' Default is 1, \emph{i.e.} 100\%
#' @param threshold A number. The TWDTW threshold, i.e. the maximum TWDTW 
#' cost for consideration. Default is \code{Inf}
#' @param normalized Use normalized TWDTW distance. Default is TRUE
#' @docType methods
#' @return object of class \code{\link[base]{data.frame}} with the best alignment 
#' for each interval
#' @examples
#' malig = mtwdtw(query.list, timeseries=template, weight = "logistic", 
#'          normalize=TRUE, query.length=23, alpha = 0.1, beta = 100)
#'          
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' best_class = classifyIntervals(x=malig, from=from, to=to, by = "6 month",
#'              overlap=.3, threshold=Inf)
#' best_class
#' 
#' 
#' @export
classfyIntervals = function(x, from, to, by, breaks=NULL, overlap=.5, threshold=Inf, normalized)
{
  .Deprecated("classifyIntervals")
  
  if (!missing(normalized))
    warning("argument normalized is deprecated and is scheduled to be removed in the next version", 
            call. = FALSE)
  
  if(is(x, "dtwSat"))
    x = getAlignments(x)
  
  if(!is(x, "data.frame"))
    stop("x is not a data.frame or dtwSat class")
  
  if( overlap < 0 & 1 < overlap )
    stop("overlap out of range, it must be a number between 0 and 1")
  
  if(is.null(breaks))
    breaks = seq(as.Date(from), as.Date(to), by=by)
  
  res = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    .bestInterval(x, start=breaks[i-1], end=breaks[i], overlap)
  }))
  
  d = res$distance
  I = d>threshold
  if(any(I)){
    res$query[I] = "unclassified"
    res$distance[I] = Inf 
  }
  res
}


.bestInterval = function(x, start, end,  overlap){
  
  I = lapply( 1:nrow(x), function(i){
    dates = seq(x$from[i], x$to[i], 1)
    dates_in = which(start <= dates & dates < end)
    r1 = length(dates_in) / as.numeric(end-start)
    if( overlap < r1 & r1 < 2-overlap )
      return(i)
    NULL
  })
  I = unlist(I)
  
  res = list()
  res$query = "unclassified"
  res$from = start
  res$to = end - 1
  res$distance = Inf
  
  if(!is.null(I)){
    i_min = which.min(x$distance[I])
    res$query = x$query[I][i_min]
    res$from = start
    res$to = end - 1
    res$distance = x$distance[I][i_min]
  }
  res = data.frame(res)
  res
}


#' @title Normalize query length 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function normalize the length of a query or 
#' list of queries
#' 
#' @param ... \link[zoo]{zoo} objects
#' @param query a list of \link[zoo]{zoo} objects
#' @param query.length An integer. Default is the length of the longest query
#' @docType methods
#' @return A \link[base]{list} of \link[zoo]{zoo} objects with normalized 
#' length 
#'
#' @examples
#' new.query.list = normalizeQuery(query = query.list, query.length = 23)
#' lapply(query.list, nrow)
#' lapply(new.query.list, nrow)
#' 
#' @export
#' 
normalizeQuery = function(..., query = list(...), query.length=NULL){
  if(is.null(query.length))
    query.length = max(unlist(lapply(query, nrow)), na.rm=TRUE)

  res = lapply(query, function(q){
    freq = as.numeric(diff(range(index(q))))/(query.length-1)
    timeline = seq(min(index(q), na.rm = TRUE), max(index(q), na.rm = TRUE), by=freq)
    na.spline(q, xout = timeline)
  })
  res
}




