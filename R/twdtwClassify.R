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
#' @inheritParams twdtwReduceTime
#'
#' @param x An object of class twdtw*. This is the target time series. 
#' Usually, it is a set of unclassified time series. 
#' 
#' @param from A character or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd".
#' 
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} object in the format "yyyy-mm-dd".
#' 
#' @param by A \code{\link[base]{character}} with the interval size, \emph{e.g.} "6 month".
#' 
#' @param breaks A vector of class \code{\link[base]{Dates}}. This replaces the arguments \code{from},
#' \code{to}, and \code{by}.
#' 
#' @param overlap A number between 0 and 1. The minimum overlapping 
#' between one match and the interval of classification. Default is 0.5, 
#' \emph{i.e.} an overlap minimum of 50\%.
#' 
#' @param patterns.labels a vector with labels of the patterns.
#' 
#' @param thresholds A numeric vector the same length as \code{patterns.labels}. 
#' The TWDTW dissimilarity thresholds, i.e. the maximum TWDTW cost for consideration 
#' in the classification. Default is \code{Inf} for all \code{patterns.labels}.
#' 
#' @param fill A character to fill the classification gaps. 
#' For signature \code{twdtwTimeSeries} the default is \code{fill="unclassified"}, 
#' for signature \code{twdtwRaster} the default is \code{fill="unclassified"}.
#' 
#' @param filepath A character. The path at which to save the raster with results. 
#' If not provided the function saves in the same directory as the input time series raster. 
#'  
#' @param ... Arguments to pass to specific methods for each twdtw* class 
#' and other arguments to pass to \code{\link[raster]{writeRaster}} and 
#' \code{\link[raster]{pbCreate}}. If \code{x} of 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}} additional arguments passed to 
#' \code{\link[dtwSat]{twdtwApply}}. 
#' 
#' @return An object of class twdtw*.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' 
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @export  
setGeneric(name = "twdtwClassify", 
          def = function(x, ...) standardGeneric("twdtwClassify"))


#' @rdname twdtwClassify 
#' @aliases twdtwClassify-data.frame 
#' @example examples/test_twdtw_raster_analysis.R 
#' @export
setMethod(f = "twdtwClassify", "data.frame",
          def = function(x, y, step.matrix=symmetric1, breaks=NULL, from=NULL, to=NULL, by=NULL, 
                         overlap=0.5,fill=length(y),alpha=-0.1,beta=50,time.window=FALSE, keep=FALSE, ...){
            twdtwReduceTime(x=x, y=y, step.matrix=step.matrix, breaks=breaks, from=from, to=to, by=by, 
                            overlap=overlap,fill=fill,alpha=alpha,beta=beta,time.window=time.window,keep=keep, ...)
          })

#' @rdname twdtwClassify 
#' @aliases twdtwClassify-list
#' @example examples/test_twdtw_raster_analysis.R 
#' @export
setMethod(f = "twdtwClassify", "list",
          def = function(x, y, step.matrix=symmetric1, breaks=NULL, from=NULL, to=NULL, by=NULL, 
                         overlap=0.5,fill=length(y),alpha=-0.1,beta=50,time.window=FALSE, keep=FALSE, ...){
            lapply(x, FUN = twdtwReduceTime, y=y, step.matrix=step.matrix, breaks=breaks, from=from, to=to, by=by, 
                   overlap=overlap,fill=fill,alpha=alpha,beta=beta,time.window=time.window,keep=keep, ...)
          })

#' @rdname twdtwClassify
#' @aliases twdtwClassify-twdtwTimeSeries 
#' @example examples/test_twdtw_raster_analysis.R 
#' @export
setMethod("twdtwClassify", "twdtwTimeSeries",
          function(x, patterns.labels=NULL, from=NULL, to=NULL, by=NULL, breaks=NULL,
                overlap=.5, thresholds=Inf, fill="unclassified", ...){
                    xm = twdtwApply(x = x, from = from, to = to, by = by, breaks = breaks, ...)
                    if(is(xm, "twdtwMatches")){
                      x = xm 
                      if(is.null(patterns.labels)) patterns.labels = labels(x@patterns)
                      if( overlap < 0 & 1 < overlap )
                        stop("overlap out of range, it must be a number between 0 and 1")
                      if(is.null(breaks))
                        if( !is.null(from) &  !is.null(to) ){
                          breaks = seq(as.Date(from), as.Date(to), by=by)    
                        } else {
                          # These automatic breaks needs to be improved 
                          y = x@patterns
                          patt_range = lapply(index(y), range)
                          patt_diff = trunc(sapply(patt_range, diff)/30)+1
                          min_range = which.min(patt_diff)
                          by = patt_diff[[min_range]]
                          cycles = c(18,12,6,4,3,2)
                          by = cycles[which.min(abs(by-cycles))]
                          from = patt_range[[min_range]][1]
                          to = from 
                          month(to) = month(to) + by
                          dates = as.Date(unlist(index(x@timeseries)))
                          year(from) = year(min(dates))
                          year(to) = year(max(dates))
                          breaks = seq(from, to, paste(by,"month"))
                        }
                      breaks = as.Date(breaks)
                      twdtwClassify.twdtwMatches(x, patterns.labels=patterns.labels, breaks=breaks, 
                                                 overlap=overlap, thresholds=thresholds, fill=fill)
                    } else {
                      new("twdtwMatches", timeseries=xm$x, patterns=xm$y, alignments=xm$aligs)
                    }
           })

#' @rdname twdtwClassify
#' @aliases twdtwClassify-twdtwTimeSeries 
#' @example examples/test_twdtw_raster_analysis.R 
#' @export
setMethod("twdtwClassify", "twdtwMatches",
          function(x, patterns.labels=NULL, from=NULL, to=NULL, by=NULL, breaks=NULL,
                   overlap=.5, thresholds=Inf, fill="unclassified"){
            if(is.null(patterns.labels)) patterns.labels = labels(x@patterns)
            if( overlap < 0 & 1 < overlap )
              stop("overlap out of range, it must be a number between 0 and 1")
            if(is.null(breaks))
              if( !is.null(from) &  !is.null(to) ){
                breaks = seq(as.Date(from), as.Date(to), by=by)    
              } else {
                # These automatic breaks needs to be improved 
                y = x@patterns
                patt_range = lapply(index(y), range)
                patt_diff = trunc(sapply(patt_range, diff)/30)+1
                min_range = which.min(patt_diff)
                by = patt_diff[[min_range]]
                cycles = c(18,12,6,4,3,2)
                by = cycles[which.min(abs(by-cycles))]
                from = patt_range[[min_range]][1]
                to = from 
                month(to) = month(to) + by
                dates = as.Date(unlist(index(x@timeseries)))
                year(from) = year(min(dates))
                year(to) = year(max(dates))
                breaks = seq(from, to, paste(by,"month"))
              }
            breaks = as.Date(breaks)
            twdtwClassify.twdtwMatches(x, patterns.labels=patterns.labels, breaks=breaks, 
                                       overlap=overlap, thresholds=thresholds, fill=fill)
          })

#' @rdname twdtwClassify
#' @aliases twdtwClassify-twdtwRaster 
#' @example examples/test_twdtw_raster_analysis.R
setMethod("twdtwClassify", "twdtwRaster",
          function(x, patterns.labels=NULL, thresholds=Inf, fill=255, filepath="", ...){
                  if(is.null(patterns.labels)) patterns.labels = coverages(x)
                  patterns.labels = patterns.labels[!patterns.labels%in%"doy"]
                  twdtwClassify.twdtwRaster(x, patterns.labels=patterns.labels, thresholds=thresholds, fill=fill, filepath=filepath, ...)
           })
           
twdtwClassify.twdtwRaster = function(x, patterns.labels, thresholds, fill, filepath, ...){
    
  if(thresholds == Inf) {
    thresholds = 9999
  }
  
  levels = c(seq_along(patterns.labels), fill)
  labels = c(patterns.labels, "unclassified")
  
  # Create output raster objects
  class_b <- brick(x@timeseries[[1]], nl = length(index(x)), values = FALSE)
  distance_b <- brick(x@timeseries[[1]], nl = length(index(x)), values = FALSE)
  class_vv <- matrix(class_b, ncol = nlayers(class_b))
  distance_vv <- matrix(distance_b, ncol = nlayers(distance_b))
  names(class_b) = paste0("date.",index(x))
  names(distance_b) = paste0("date.",index(x))

  filepath <- trim(filepath)
  filename <- NULL
  if (filepath != "") {
    dir.create(path = filepath, showWarnings = TRUE, recursive = TRUE)
    filename <- paste0(filepath, "/", c("Class", "Distance"), ".grd")
    names(filename) <- c("Class", "Distance")
  } else if (!canProcessInMemory(class_b, n = length(x@timeseries) + 2)) {
    filename <- c(rasterTmpFile("Class"), rasterTmpFile("Distance"))
  }

  if (!is.null(filename)) {
    class_b <- writeStart(class_b, filename = filename[1], ...)
    distance_b <- writeStart(distance_b, filename = filename[2], ...)
  }

  bs <- blockSize(x@timeseries[[1]])
  bs$array_rows <- cumsum(c(1, bs$nrows * class_b@ncols))
  pb <- pbCreate(bs$n, ...)

  for(k in 1:bs$n){
    
    v <- lapply(x@timeseries[patterns.labels], getValues, row = bs$row[k], nrows = bs$nrows[k])
    rows <- seq(from = bs$array_rows[k], by = 1, length.out = bs$nrows[k]*class_b@ncols)
    
    for(i in seq_along(index(x))) {
      
      r <- sapply(v, function(vv) vv[, i])
      d <- apply(r, 1, min)
      dc <- apply(r, 1, which.min)
      dc[which(d[]>=thresholds)] = fill
      class_vv[rows, i] <- dc
      distance_vv[rows, i] <- d
    }
    
    if (!is.null(filename)) {
      writeValues(class_b, class_vv[rows, ], bs$row[k])
      writeValues(distance_b, distance_vv[rows, ], bs$row[k])
    } 
    
    pbStep(pb, k)
    
  }
  
  if (!is.null(filename)) {
    class_b <- writeStop(class_b)
    distance_b <- writeStop(distance_b)
  } else {
    class_b <- setValues(class_b, values = class_vv)
    distance_b <- setValues(distance_b, values = distance_vv)
  }
  
  pbClose(pb)
  
  twdtwRaster(Class = class_b, Distance = distance_b, ..., timeline = index(x), 
              labels = labels, levels = levels, filepath = filepath)
                
}

twdtwClassify.twdtwMatches = function(x, patterns.labels, breaks, overlap, thresholds, fill){
    levels = as.character(patterns.labels)
    names(levels) = levels
    m = length(levels)
    n = length(breaks)-1
    aligs = lapply(as.list(x), FUN=.bestIntervals, m=m, n=n, levels=levels, breaks=breaks, overlap=overlap)
    twdtwMatches(x@timeseries, patterns=x@patterns, alignments=aligs)
}

.bestIntervals = function(x, m, n, levels, breaks, overlap)
{
  best_matches = .bestmatches(x, m, n, levels, breaks, overlap)$IM
  IL = best_matches[,1]
  I = unique(best_matches[,1])
  I = I[I>0]
  names(I) = levels[I]
  aligs = lapply(levels, initAlignments)
  aligs[names(I)] = lapply(I, function(i) subset(x, timeseries.labels = 1, patterns.labels = i, k = best_matches[IL==i,3])@alignments[[1]][[1]] )
  new("twdtwMatches", x@timeseries, x@patterns, alignments=list(aligs))
}
