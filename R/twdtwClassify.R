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
#' @name twdtwClassify
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
#' @param thresholds A numeric vector the same length as \code{patterns.labels}. 
#' The TWDTW dissimilarity thresholds, i.e. the maximum TWDTW cost for consideration 
#' in the classification. Default is \code{Inf} for all \code{patterns.labels}.
#' 
#' @param fill a character or value to fill the classification gaps. 
#' For signature \code{twdtwTimeSeries} the default is \code{fill="unclassified"}, and 
#' for signature \code{twdtwRaster} the default is \code{fill="unclassified"}.
#'  
#' @param ... arguments to pass to specifique methods for each twdtw* signature 
#' and other arguments to pass to \code{\link[raster]{writeRaster}}. 
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
setGeneric(name = "twdtwClassify", 
          def = function(x, ...) standardGeneric("twdtwClassify"))

#' @rdname twdtwClassify
#' @aliases twdtwClassify-twdtwTimeSeries 
#' @examples
#' # Classifying time series based on TWDTW results 
#' ts = twdtwTimeSeries(example_ts.list)
#' patt = twdtwTimeSeries(patterns.list)
#' log_fun = logisticWeight(-0.1, 100)
#' time_intervals = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), by="6 month")
#' mat = twdtwApply(x=ts, y=patt, weight.fun=log_fun)
#' best_mat = twdtwClassify(mat, breaks=time_intervals, overlap=0.5)
#' 
#' \dontrun{
#' require(parallel)
#' best_mat = mclapply(mat, mc.cores=4, FUN=twdtwClassify, breaks=time_intervals, overlap=0.5) 
#' best_mat = twdtwMatches(alignments=best_mat)
#' }
#' @export
setMethod("twdtwClassify", "twdtwMatches",
          function(x, patterns.labels=NULL, from=NULL, to=NULL, by=NULL, breaks=NULL,
                overlap=.5, thresholds=Inf, fill="unclassified"){
                    if(is.null(patterns.labels)) patterns.labels = labels(x@patterns)
                    if( overlap < 0 & 1 < overlap )
                      stop("overlap out of range, it must be a number between 0 and 1")
                    if(is.null(breaks))
                      breaks = seq(as.Date(from), as.Date(to), by=by)
                    breaks = as.Date(breaks)
                  twdtwClassify.twdtwMatches(x, patterns.labels=patterns.labels, breaks=breaks, 
                            overlap=overlap, thresholds=thresholds, fill=fill)
           })

#' @rdname twdtwClassify
#' @aliases twdtwClassify-twdtwRaster 
#' @examples
#' \dontrun{
#' # Run TWDTW analysis for raster time series 
#' load(system.file("lucc_MT/temporal_patterns.RData", package="dtwSat"))
#' patt = twdtwTimeSeries(temporal_patterns)
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' 
#' time_interval = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), 
#'                     by="6 month")
#' log_fun = weight.fun=logisticWeight(-0.1,100)
#' 
#' r_twdtw = twdtwApply(x=rts, y=patt, weight.fun=log_fun, breaks=time_interval, 
#'           filepath="~/test_twdtw", overwrite=TRUE, format="GTiff", mc.cores=3)
#' 
#' lucc = twdtwClassify(r_twdtw, format="GTiff")
#' 
#' }
setMethod("twdtwClassify", "twdtwRaster",
          function(x, patterns.labels=NULL, thresholds=Inf, fill=255, ...){
                  if(is.null(patterns.labels)) patterns.labels = names(x)[-1]
                  twdtwClassify.twdtwRaster(x, patterns.labels=patterns.labels, thresholds=thresholds, fill=fill, ...)
           })
           
twdtwClassify.twdtwRaster = function(x, patterns.labels, thresholds, fill, ...){
    
    if(thresholds==Inf) thresholds = 9999
        
    out = lapply(seq_along(index(x)), function(i) {
      r = lapply(as.list(x)[patterns.labels], raster, layer=i)
      b = brick(r)
      mb = min(b)
      res = which.min(b)
      res[which(mb[]>=thresholds)] = fill
      list(class=res, distance=mb)
    })
    
    class_b = do.call("brick", lapply(out, function(x) x$class))
    distance_b = do.call("brick", lapply(out, function(x) x$distance)) 
    names(class_b) = paste0("date.",index(x)) 
    names(distance_b) = paste0("date.",index(x)) 
    
    labels = factor( values, levels = values, labels = c(patterns.labels, "unclassified") )
    
    if(is.null(filepath)) filepath = filename(x[[2]])
    twdtwRaster(Class=class_b, Distance=distance_b, ..., timeline=index(x), 
                labels=labels, filepath=dirname(filepath))
                
}

twdtwClassify.twdtwMatches = function(x, patterns.labels, breaks, overlap, thresholds, fill){
    aligs = lapply(as.list(x), FUN = classifyIntervals, patterns.labels, breaks, overlap, thresholds, fill)
    twdtwMatches(x@timeseries, patterns=x@patterns, alignments=aligs)
}

classifyIntervals = function(x, patterns.labels, breaks, overlap, thresholds, fill)
{
  
  dist_table = x[[1, patterns.labels]]
  labels = as.character(patterns.labels)
  names(labels) = labels
  
  best_match = do.call("rbind", lapply(seq_along(breaks)[-1], function(i){
    from = breaks[i-1]
    to = breaks[i]
    data.frame(from, to, K= .bestMatches(x=dist_table, start=from, end=to, overlap))
  }))
  best_match = best_match[!is.na(best_match$K),]
  
  if(nrow(best_match)<1) {
      res = lapply(names(labels), function(p){
            alignments = list()
            alignments$label = numeric(0)
            alignments$from = numeric(0)
            alignments$to = numeric(0)
            alignments$distance = numeric(0)
            alignments$K = 0
            alignments$matching = list()
            alignments$internals = getInternals(x, 1, p)[[1]][[p]]$internals
            alignments
      })
  } else {
  
      K = best_match$K
      
      best_match = best_match[ dist_table$distance[best_match$K]<=thresholds , ]
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
      res = res[sapply(res, function(r) r$K>0)]
  }
  names(res) = names(labels)
  res
}

.bestMatches = function(x, start, end, overlap){
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

.lowestDistances = function(x, start, end, overlap){
  J = (x$from <= end & x$to >= start)
  x = x[J,]
  x$from[x$from < start] = start
  x$to[end < x$to] = end
  # Check for minimum overlap 
  r1 = as.numeric(x$to - x$from) / as.numeric(end-start)
  I = overlap < r1 & r1 < 2-overlap
  J[which(J)] = I
  if(!any(I)) return(9999)
  # Sellect the lowest TWDTW distance 
  min(x$distance[I])
}