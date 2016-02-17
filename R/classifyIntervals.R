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
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @title Classify time intervals
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the best match within each 
#' predefined period based on the TWDTW dissimilarity measure.
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object.
#' 
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd".
#' 
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd".
#' 
#' @param by A \code{\link[base]{character}} with the intevals size, \emph{e.g.} "6 month".
#' 
#' @param breaks A vector of class \code{\link[base]{Dates}}. This replaces the arguments \code{from},
#' \code{to}, and \code{by}.
#' 
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between one match and the interval of classification. Default is 0.5, 
#' \emph{i.e.} an overlap minimum of 50\%.
#' 
#' @param threshold A number. The TWDTW dissimilarity threshold, i.e. the maximum TWDTW 
#' cost for consideration in the classification. Default is \code{Inf}.
#' 
#' @param simplify A logical. TRUE returns a vector with the best pattern for 
#' each interval. FALSE returns a \code{\link[dtwSat]{twdtw-class}} object. 
#' Default is FALSE.
#'  
#' @param y A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' considers all patterns in the classification. 
#' 
#' @docType methods
#' 
#' @return A \code{\link[dtwSat]{twdtw-class}} with the best match 
#' for each classification period or a vector with the best class 
#' for each period.
#' 
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' 
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun)
#'          
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' by = "6 month"
#' 
#' # All classes
#' class_1 = classifyIntervals(x=matches, from=from, to=to, by = by,
#'              overlap=.5)
#' 
#' getAlignments(class_1)
#' plotAlignments(class_1)
#' plotClassification(class_1)
#' 
#' # Cotton and Maize 
#' class_2 = classifyIntervals(x=matches, y=c("Cotton","Maize"), 
#'                  from=from, to=to, by = by, overlap=.5)
#' 
#' getAlignments(class_2)
#' plotAlignments(class_2)
#' plotClassification(class_2)
#' 
#' # Simplify Cotton and Maize 
#' classifyIntervals(x=matches, from=from, to=to, by = by,
#'              overlap=.5, threshold=Inf, simplify=TRUE)
#' 
#'         
#' @export
classifyIntervals = function(x, y, from=NULL, to=NULL, by=NULL, breaks=NULL,
                             overlap=.5, threshold=Inf, simplify=FALSE)
{
  
  y = getPatternNames(x, y)

  if(!is(x, "twdtw"))
    stop("aligs is not twdtw-class")
  
  if( overlap < 0 & 1 < overlap )
    stop("overlap out of range, it must be a number between 0 and 1")
  
  if(is.null(breaks))
    breaks = seq(as.Date(from), as.Date(to), by=by)
  
  breaks = as.Date(breaks)
  
  aligs = getAlignments(x, y)
  I = sapply(seq_along(breaks)[-1], function(i){
    .bestMatche(aligs, start=breaks[i-1], end=breaks[i], overlap)
  })
  
  d = aligs$distance[I]
  I[ d>threshold ] = NA
  res = aligs$pattern[I]
  names(res) = paste0("Period", seq_along(breaks[-1]))
  
  if(!simplify){
    J = !is.na(I)
    aux_matches = getMatches(x, y)
    aux_internals = getInternals(x, y)
    n_matches = sapply(aux_matches, length)
    P = unlist(sapply(sapply(aux_matches, length), seq, from=1))
    S = P[seq_along(P)[I[J]]]
    best_matches = lapply(y, function(p){
      L = S[aligs$pattern[I[J]]==p]
      alignments = list()
      alignments$pattern = res[J][res[J]==p]
      alignments$from = breaks[-length(breaks)][J][aligs$pattern[I[J]]==p]
      alignments$to = breaks[-1][J][aligs$pattern[I[J]]==p]
      alignments$distance = d[J][res[J]==p]
      alignments$K = length(L)
      if(n_matches[p]>0){
        alignments$matching = aux_matches[[p]][L]
        alignments$internals = aux_matches[[p]][L]
      }
      alignments
    })
    res = new("twdtw", call=match.call(), x = getTimeSeries(x), patterns = getPatterns(x, y), alignments = best_matches)
  }
  res 
}


.bestMatche = function(x, start, end, overlap){
  J = (x$from <= end & x$to >= start)
  x = x[J,]
  x$from[x$from < start] = start
  x$to[end < x$to] = end
  # Check for minimum overlap 
  r1 = as.numeric(x$to - x$from) / as.numeric(end-start)
  I = overlap < r1 & r1 < 2-overlap
  J[which(J)] = I
  if(!any(I)) return(NA)
  # Sellect the lowest TWDTW distance 
  res = which(J)[which.min(x$distance[I])]
  res
}

