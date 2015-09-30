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
#### UTILITY FUNCTIONS



#' @title Wavelet filter
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a smoothing algorithm to
#' the time series. It computes a discreat wavelet 
#' smoothing for each dimension in the imput time series.
#' 
#' @param x A \code{\link[zoo]{zoo}} object with the time series
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
#' @docType methods
#' @return object of class \code{\link[zoo]{zoo}} 
#' 
#' @seealso \link[waveslim]{mra}
#' 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' plot(template$evi, ylab="EVI", xlab="Time")
#' lines(sy$evi, col="red")
#' 
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
#' @return vector of class \code{\link[base]{Date}} 
#' 
#' @seealso \link[dtwSat]{getDatesFromDOY} and \link[dtwSat]{getModisTimeIndex}
#' 
#' @examples
#' dates = createTimeSequence()
#' dates
#' 
#' @export
createTimeSequence = function(year=2000:format(Sys.time(), "%Y"), frequency=16){
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
#' @return object of class \code{\link[base]{Date}} 
#' 
#' @seealso \link[dtwSat]{createTimeSequence} and \link[dtwSat]{getModisTimeIndex}
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
#' @param date A \code{\link[base]{Date}} object
#' @param frequency An integer with the frequency in days. Default is 16 days
#' 
#' @docType methods 
#' 
#' @return An integer 
#' 
#' @seealso \link[dtwSat]{createTimeSequence} and \link[dtwSat]{getDatesFromDOY}
#' 
#' @examples
#' i = getModisTimeIndex(date=as.Date("2000-01-01"))
#' i
#'
#' @export
getModisTimeIndex = function(date, frequency=16){
  dates = createTimeSequence(frequency=frequency)
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
#' @param from A \code{\link[base]{Date}} object
#' @param by A character with the intevals size, \emph{e.g.} ''6 month''
#' @param to A \code{\link[base]{Date}} object
#' @param breaks A vector of class \code{\link[base]{Date}}
#' @param normalized Use normalized TWDTW distance. Default is TRUE
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between the one alignment and the interval of classification. 
#' Default is 1, \emph{i.e.} 100\%
#' @param threshold A number. The TWDTW threshold, i.e. the maximum TWDTW 
#' cost for consideration. Default is \code{Inf}
#' @docType methods
#' @return object of class \code{\link[base]{data.frame}} with the best alignment 
#' for each interval
#' @examples
#' malig = mtwdtw(query.list, template, weight = "logistic", 
#'          alpha = 0.1, beta = 100)
#'          
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' best_class = classfyIntervals(x=malig, from=from, to=to, by = "6 month",
#'              normalized=TRUE, overlap=.7, threshold=Inf)
#' best_class
#' 
#' 
#' @export
classfyIntervals = function(x, from, to, by, breaks=NULL,
                            normalized=TRUE, overlap=1, threshold=Inf)
{
  if(is(x, "dtwSat"))
    x = getAlignments(x)

  if(!is(x, "data.frame"))
    stop("x is not a data.frame or dtwSat class")
  
  if( overlap < 0 & 1 < overlap )
    stop("overlap out of range, it must be a number between 0 and 1")

  if(is.null(breaks))
    breaks = seq(from, to, by=by)
  
  res = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    .bestInterval(x, start=breaks[i-1], end=breaks[i], normalized, overlap)
  }))
  
  d = res$distance
  if(normalized)
    d = res$normalizedDistance
  res[d<=threshold,]
}

.bestInterval = function(x, start, end, normalized, overlap){
  # Find alignments within the interval
  in_interval = lapply( 1:nrow(x), function(i){
    dates = seq(x$from[i], x$to[i], 1)
    dates_in = which(start <= dates & dates < end)
    if( length(dates_in) / length(dates) < overlap )
      return(NA)
    i
  })

  # Find lowest cost alignment
  in_interval = unlist(in_interval)
  i_min = which.min(x$distance[in_interval])
  if(normalized)
    i_min = which.min(x$normalizedDistance[in_interval])

  if(length(i_min)!=1)
    return(NULL)
  i_min = in_interval[i_min]
  res = list()
  res$query = x$query[i_min]
  res$from = start
  res$to = end - 1
  res$distance = x$distance[i_min]
  res$normalizedDistance = x$normalizedDistance[i_min]
  res = data.frame(res)
  res
}


#### CREATE DATASET FOR CLASSIFICATION
#### CREATE A SMALL CLASSIFIED MAP
#### FROM PAPER IEEE (Small cut of Proto dos gauchos)








