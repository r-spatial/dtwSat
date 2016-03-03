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


#' @include methods.R
#' @title Classify time series 
#' @name classifyTimeSeries
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function classifies the intervals of a time series 
#' based on the TWDTW results. 
#' 
#' @inheritParams get
#'
#' @param x an object of class twdtw*. This is the target time series. 
#' Usually, it is a set of unclassified time series. 
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
#' @param fill a character or value to fill the classification gaps. Default is \code{fill="unclassified"}.
#' 
#' @param simplify logical. TRUE returns a list of labels with the best pattern for 
#' each period of time defined by \code{breaks}. FALSE returns an object of class twdtwMatches 
#' with the best matches. Default is FALSE. 
#' 
#'  
#' @param ... arguments to pass to specifique methods for each twdtw* signature 
#' and other arguments for parallel processing to pass
#'
#' @return An object of class twdtw*.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' 
#' @export  
setGeneric(name = "classifyTimeSeries", 
          def = function(x, patterns.labels=NULL, from=NULL, to=NULL, by=NULL, breaks=NULL,
                overlap=.5, threshold=Inf, fill="unclassified", simplify=FALSE, ...) standardGeneric("classifyTimeSeries"))


#' @rdname classifyTimeSeries
#' @aliases classifyTimeSeries-twdtwTimeSeries 
#' @examples
#' # Applying TWDTW analysis to objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(example_ts)
#' patt = twdtwTimeSeries(patterns.list)
#' mat = twdtwApply(x=ts, y=patt, weight.fun=logisticWeight(-0.1, 100), keep=TRUE)
#' best_mat = classifyTimeSeries(mat, from="2009-09-01", to="2013-09-01", by="6 month", overlap=0.5) 
#'
#' @export
setMethod("classifyTimeSeries", "twdtwMatches",
          function(x, patterns.labels, from, to, by, breaks, overlap, threshold, fill, simplify, ...){
                    if(is.null(patterns.labels)) patterns.labels = labels(x@patterns)
                    if( overlap < 0 & 1 < overlap )
                      stop("overlap out of range, it must be a number between 0 and 1")
                    if(is.null(breaks))
                      breaks = seq(as.Date(from), as.Date(to), by=by)
                    breaks = as.Date(breaks)
                  classifyTimeSeries.twdtwMatches(x, patterns.labels=patterns.labels, breaks=breaks, 
                                                  overlap=overlap, threshold=threshold, fill=fill, simplify=simplify, ...)
           })
         
classifyTimeSeries.twdtwMatches = function(x, patterns.labels, breaks, overlap, threshold, fill, simplify, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 1L), mc.cleanup = TRUE){
    if(mc.cores>1){
      aligs = mclapply(as.list(x), mc.preschedule=mc.preschedule, mc.set.seed=mc.set.seed, mc.silent=mc.silent, 
                                 mc.cores=mc.cores, mc.cleanup=mc.cleanup, FUN = classifyIntervals, patterns.labels, breaks, overlap, threshold, fill, simplify)
    } else {
      aligs = lapply(as.list(x), FUN = classifyIntervals, patterns.labels, breaks, overlap, threshold, fill, simplify)
    }
    if(simplify) return(aligs)
    twdtwMatches(x@timeseries, patterns=x@patterns, alignments=aligs)
}


classifyIntervals = function(x, patterns.labels, breaks, overlap, threshold, fill, simplify){
  
  dist_table = x[[1, patterns.labels]]
  labels = as.character(patterns.labels)
  names(labels) = labels
  
  best_match = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    from = breaks[i-1]
    to = breaks[i]
    data.frame(from, to, K= .bestMatche(x=dist_table, start=from, end=to, overlap))
  }))
  
  K = best_match$K
    
  if(simplify) {
    res = as.character(dist_table$label[K])
    res[dist_table$distance[best_match$K]>threshold] = fill
    names(res) = paste0("Period",seq_along(res))
  } else {
    best_match = best_match[ dist_table$distance[best_match$K]<=threshold , ]
    K = best_match$K
    P = lapply(labels, function(p) K[dist_table$label[K]==p])
    from = lapply(labels, function(p) best_match$from[dist_table$label[K]==p])
    to = lapply(labels, function(p) best_match$to[dist_table$label[K]==p])
    J = unlist(lapply(x@alignments[[1]], seq_along))
    res = lapply(names(P), function(p){
        alignments = list()
        alignments$label = p
        alignments$from = from[[p]]
        alignments$to = to[[p]]
        alignments$distance = dist_table$distance[P[[p]]]
        alignments$K = length(P[[p]])
        alignments$matching = lapply(P[[p]], function(k) { 
            matching = getMatches(x, 1, p)[[1]][[p]]$matching
            if(length(matching)<J[k]) return(NULL)
            matching[[J[k]]]
        })
        alignments$matching = alignments$matching[!sapply(alignments$matching, is.null)]
        alignments$internals = getInternals(x, 1, p)[[1]][[p]]$internals
        alignments
    })
    names(res) = names(P)
    res = res[sapply(res, function(r) r$K>0)]
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

