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
#   R Package dtwSat - 2016-16-01                             #
#                                                             #
###############################################################


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
#' @param simplify return only the best pattenr for each interval. Default is FALSE 
#' @param levels A character or numeric vector. The categories for classification
#' @param labels A character or numeric vector. The labels for each category
#' @param Unclassified A numeric to fill gaps. Default is 255
#' @param ... other argument passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' 
#' @docType methods
#' @return A \code{\link[base]{data.frame}} with the best alignment 
#' for each interval
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
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
#' # Simplify Cotton and Maize 
#' classifyIntervals(x=alig, from=from, to=to, by = by, simplify= TRUE,
#'              overlap=.3, threshold=Inf, p.names=c("Cotton","Maize"))
#'              
#' @export
classifyIntervals = function(x, breaks=NULL, from=NULL, to=NULL, by=NULL,
                             overlap=.3, threshold=Inf, 
                             simplify=FALSE,
                             levels=NULL,
                             labels=NULL,
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
  
  if(is.null(levels))
    levels = c(seq_along(unique(x$pattern)), Unclassified)
  
  if(is.null(labels))
    labels = c(unique(x$pattern), "Unclassified")
  
  if(!length(levels)==length(labels))
    stop("levels and labels are not the same length")
  
  I = match(res$pattern, labels)
  res$level = levels[I]
  names(res$level) = labels[I]
  
  if(simplify)
    return(res$level)
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

