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

#' @title Create raster time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a list of raster time series 
#' 
#' @param ... \code{\link[raster]{RasterBrick-class}} or \code{\link[raster]{RasterStack-class}}
#' @param x A list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}}
#' Each layer of the Raster* object is a time step
#' @param timeline A vector of \code{\link[base]{Dates}}
#' It must have the length of the layers in the \code{\link[raster]{RasterBrick-class}} object 
#' @param doy optional. A \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} with the same extent as the 
#' objects in \code{x}
#' 
#' @details a list of \code{\link[raster]{RasterBrick-class}} or 
#' \code{\link[raster]{RasterStack-class}} objects
#' 
#' @docType methods
#' @return A \code{\link[dtwSat]{twdtw-class}} object
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @examples
#' 
#' ####
#' 
#' @export
createRasterTS = function(... , x = list(...), timeline, doy){
  
  timeline = as.Date(timeline)
  if(!exists("doy")) {
    array_data = rep(as.numeric(format(timeline, "%j")), each=ncell(x[[1]]))
    doy = array(array_data, dim = dim(x[[1]]))
  }
  x = c(x, doy=doy)
  
  # Check timeline 
  nl = lapply(x, nlayers)
  if(any(nl!=length(timeline)))
    stop("raster objects do not have the same length as the timeline")
  
  fun = function(x){
    names(x) = paste0("date.",timeline)
    x
  }
  
  res = lapply(x, FUN=fun)
  res
}

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
#' @param ... other arguments to pass to the function \code{\link[waveslim]{mra}} in the 
#' packege \pkg{waveslim}
#' 
#' @docType methods
#' @return A \code{\link[zoo]{zoo}} object
#' 
#' @seealso \link[waveslim]{mra}
#' 
#' @examples
#' ## Wavelet filter
#' sy = waveletSmoothing(x=template, frequency=16, wf = "la8", J=1, 
#'      boundary = "periodic")
#' 
#' ## Plot raw EVI and filtered EVI
#' evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
#' gp = autoplot(evi, facets = NULL)
#' gp
#'     
#' ## Plot all bands
#' gp = autoplot(sy, facets = NULL)
#' gp
#' 
#' @export
waveletSmoothing = function(x, timeline=NULL, frequency=NULL, 
                            wf = "la8", J=1, boundary = "periodic", ...)
{

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
  y = zoo(order.by = timeline)
  y = merge(x, y)
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
#' @return A vector of \code{\link[base]{Dates}} 
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
#' @return A \code{\link[base]{Dates}} object
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
  res = as.Date(paste(as.numeric(year), as.numeric(doy)), format="%Y %j", origin="1970-01-01")
  # Correct leap years 
  I = which(diff(res)<0)
  if(length(I)>0) res[I+1] = as.Date(paste0(as.numeric(format(res[I+1], "%Y"))+1,format(res[I+1],"-%m-%d")))
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
#' @return An \code{\link[base]{integer}} object 
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
#' @param x A \code{\link[dtwSat]{twdtw-class}} object or 
#' a \code{\link[base]{data.frame}} such as retrieved by \code{\link[dtwSat]{getAlignments}} 
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd"
#' @param by A \code{\link[base]{character}} with the intevals size, \emph{e.g.} ''6 month''
#' @param breaks A vector of class \code{\link[base]{Dates}}
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between the one alignment and the interval of classification. 
#' Default is 1, \emph{i.e.} 100\%
#' @param threshold A number. The TWDTW threshold, i.e. the maximum TWDTW 
#' cost for consideration. Default is \code{Inf}
#' @param pattern.only return only the best pattenr for each interval. Default is FALSE 
#' @param pattern.levels A character or numeric vector. The categories for classification
#' @param pattern.labels A character or numeric vector. The labels for each category
#' @param Unclassified A numeric to fill gaps. Default is 255
#' @param ... other argument passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' 
#' @docType methods
#' @return A \code{\link[base]{data.frame}} with the best alignment 
#' for each interval
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' 
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23)
#'          
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' by = "6 month"
#' 
#' # All classes
#' classifyIntervals(x=alig, from=from, to=to, by = by,
#'              overlap=.3, threshold=Inf)
#' 
#' # Cotton and Maize 
#' classifyIntervals(x=alig, from=from, to=to, by = by,
#'              overlap=.3, threshold=Inf, p.names=c("Cotton","Maize"))
#' 
#' 
#' @export
classifyIntervals = function(x, breaks=NULL, from=NULL, to=NULL, by=NULL,
                             overlap=.3, threshold=Inf, 
                             pattern.only=FALSE,
                             pattern.levels=NULL,
                             pattern.labels=NULL,
                             Unclassified=255,
                             ...)
{
  
  p.names = getPatternNames(x, ...)
  
  if(is(x, "twdtw"))
    x = getAlignments(x, p.names)

  if(!is(x, "data.frame"))
    stop("x is not a data.frame or twdtw-class")
  
  if( overlap < 0 & 1 < overlap )
    stop("overlap out of range, it must be a number between 0 and 1")

  if(is.null(breaks))
    breaks = seq(as.Date(from), as.Date(to), by=by)
  
  if(!is(breaks,"Dates"))
    breaks = as.Date(breaks)
  
  res = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    .bestInterval(x, start=breaks[i-1], end=breaks[i], overlap, Unclassified)
  }))
  
  d = res$distance
  I = d>threshold
  if(any(I)){
    res$pattern[I] = if(is(x$pattern, "character")){"Unclassified"}else{Unclassified}
    res$distance[I] = Inf 
  }
  
  if(!is.null(pattern.levels)&&!is.null(pattern.labels)){
    if(!length(pattern.levels)==length(pattern.labels))
      stop("pattern.levels and pattern.labels are not the same length")
    I = match(res$pattern, pattern.levels)
    res$pattern = pattern.labels[I]
    names(res$pattern) = pattern.levels[I]
  }
  
  if(pattern.only)
    return(res$pattern)
  res
}


.bestInterval = function(x, start, end,  overlap, Unclassified){
  
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
  res$pattern = if(is(x$pattern, "character")){"Unclassified"}else{Unclassified}
  res$from = start
  res$to = end - 1
  res$distance = Inf
  
  if(!is.null(I)){
    i_min = which.min(x$distance[I])
    res$pattern = x$pattern[I][i_min]
    res$from = start
    res$to = end - 1
    res$distance = x$distance[I][i_min]
  }
  res = data.frame(res, stringsAsFactors = FALSE)
  res
}


#' @title Normalize patterns length 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function normalizes the length of a pattern or 
#' list of patterns
#' 
#' @param ... \link[zoo]{zoo} objects
#' @param patterns a list of \link[zoo]{zoo} objects
#' @param patterns.length An integer. Default is the length of the longest pattern
#' @docType methods
#' @return A \link[base]{list} of \link[zoo]{zoo} objects with normalized 
#' length 
#'
#' @examples
#' new.patterns.list = normalizePatterns(patterns = patterns.list, patterns.length = 23)
#' sapply(patterns.list, nrow)
#' sapply(new.patterns.list, nrow)
#' 
#' @export
#' 
normalizePatterns = function(..., patterns = list(...), patterns.length=NULL){
  if(is.null(patterns.length))
    patterns.length = max(unlist(lapply(patterns, nrow)), na.rm=TRUE)

  res = lapply(patterns, function(q){
    freq = as.numeric(diff(range(index(q))))/(patterns.length-1)
    timeline = seq(min(index(q), na.rm = TRUE), max(index(q), na.rm = TRUE), by=freq)
    res = zoo(data.frame(na.spline(q, xout = timeline)), timeline)
    names(res) = names(q)
    res
  })
  res
}

#' @title Shift dates of zoo object 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function shifts the dates of the time series to a base year 
#' 
#' 
#' @param x a \link[zoo]{zoo} object 
#' @param year An integer
#' 
#' @docType methods
#' @return a \link[zoo]{zoo} object shifted to the given year 
#'
#' @examples
#' dates = seq(from = as.Date("2014-09-01"), to = as.Date("2015-03-01"), by = "15 day")
#' x = zoo(seq_along(dates), dates)
#' y = shiftDate(x, year=2005)
#' index(x)
#' index(y)
#' 
#' @export
#' 
shiftDate = function(x, year){
  doy = as.numeric(format(index(x), "%j"))
  d = .shiftDate(year, doy)
  index(x) = d
  x
}

.shiftDate = function(year, doy){
  dates = getDatesFromDOY(year, doy)
  diffdate = diff(dates)
  I = which(diffdate<0)
  if(length(I)<1) return(dates)
  dates = c(.shiftDate(year-1, doy[1:(I[1])]), dates[(I[1]+1):length(doy)])
  dates
}


# Match and set a list of arguments to a function 
.setFunArgs = function(fun, ..., args = list(...)){ 
  base_formals = formals(fun)
  base_formals_names = names(base_formals)
  given_formals = args[names(args) %in% base_formals_names]
  missing_formals_names = setdiff(base_formals_names, names(args))
  new_formals = c(base_formals[missing_formals_names], given_formals)
  new_formals = new_formals[base_formals_names]
  formals(fun) = new_formals
  fun
}



# Crop raster time series. Returns a 3D array 
.cropTimeSeries = function(x, r1, r2){
  if(is(x, "RasterBrick")){
    y = extent(x, r1, r2)
    x = crop(x, y)
  } else {
    x = x[r1:r2,,]
  }
  alply(x, c(1,2), as.numeric)
}

# Build zoo time series  
.bulidZooFromTS = function(p, x, timeline, bands){
  # Get time series for each band 
  datasets = lapply(x, function(x) x[[p]])
  datasets$doy = getDatesFromDOY(doy=datasets$doy, year=format(timeline, "%Y"))
  
  # Remove invalid values 
  k = unlist(lapply(datasets[bands], function(x){
    which(x<0|is.na(x))
  }))
  k = c(k, which(duplicated(datasets$doy)))
  if(length(k)>0) datasets = lapply(datasets, function(x) x[-k] )
  
  # Build multi-band zoo object 
  zoo(data.frame(datasets[bands]), order.by = datasets$doy)
}
