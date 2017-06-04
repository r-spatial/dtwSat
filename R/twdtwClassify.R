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
#' @param filepath A character. The path to save the raster with results. If not informed the 
#' function saves in the same directory as the input time series raster. 
#'  
#' @param ... arguments to pass to specifique methods for each twdtw* class 
#' and other arguments to pass to \code{\link[raster]{writeRaster}} and 
#' \code{\link[raster]{pbCreate}}. 
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
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' log_fun = logisticWeight(-0.1, 100)
#' time_intervals = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), by="6 month")
#' mat = twdtwApply(x=ts, y=patt, weight.fun=log_fun, keep=TRUE)
#' best_mat = twdtwClassify(x=mat, breaks=time_intervals, overlap=0.5)
#' plot(x=best_mat, type="classification")
#' 
#' \dontrun{
#' require(parallel)
#' best_mat = mclapply(as.list(mat), mc.cores=2, FUN=twdtwClassify, breaks=time_intervals, overlap=0.5)
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
                    if( !is.null(from) &  !is.null(to) ){
                      breaks = seq(as.Date(from), as.Date(to), by=by)    
                    } else {
                      # This automatic breaks needs to be improved 
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
#' @examples
#' \dontrun{
#' # Create raster time series
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
#' # Read fiels samples 
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", 
#'                 package="dtwSat"), what = "character")
#' 
#' # Split samples for training (10%) and validation (90%) using stratified sampling 
#' library(caret) 
#' set.seed(1)
#' I = unlist(createDataPartition(field_samples$label, p = 0.1))
#' training_samples = field_samples[I,]
#' validation_samples = field_samples[-I,]
#' 
#' # Get time series form raster
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' validation_ts = getTimeSeries(rts, y = validation_samples, proj4string = proj_str)
#' 
#' # Create temporal patterns 
#' temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))
#' 
#' # Set TWDTW weight function 
#' log_fun = weight.fun=logisticWeight(-0.1, 50)
#' 
#' # Run serial TWDTW analysis 
#' r_twdtw <- twdtwApply(x = rts, y = temporal_patterns, weight.fun = log_fun)
#'                                 
#' # Run parallel TWDTW analysis
#' beginCluster()
#' r_twdtw <- twdtwApplyParallel(x = rts, y = temporal_patterns, weight.fun = log_fun)
#' endCluster()
#' 
#' # Classify raster based on the TWDTW analysis 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' plot(r_lucc)
#' 
#' 
#' }
setMethod("twdtwClassify", "twdtwRaster",
          function(x, patterns.labels=NULL, thresholds=Inf, fill=255, filepath="", ...){
                  if(is.null(patterns.labels)) patterns.labels = coverages(x)
                  patterns.labels = patterns.labels[!patterns.labels%in%"doy"]
                  # if(missing(filepath)) filepath = if(fromDisk(x[[2]])){dirname(filename(x[[2]]))}else{NULL}
                  twdtwClassify.twdtwRaster(x, patterns.labels=patterns.labels, thresholds=thresholds, fill=fill, filepath=filepath, ...)
           })
           
twdtwClassify.twdtwRaster = function(x, patterns.labels, thresholds, fill, filepath, ...){
    
  if(thresholds == Inf) {
    thresholds = 9999
  }
  
  # # Create output raster objects 
  # class_b <- brick(x@timeseries[[1]], nl = length(index(x)), values = FALSE)
  # distance_b <- brick(x@timeseries[[1]], nl = length(index(x)), values = FALSE)
  # names(class_b) = paste0("date.",index(x)) 
  # names(distance_b) = paste0("date.",index(x)) 
  # 
  # filepath <- trim(filepath)
  # filename <- NULL
  # if (filepath != "") {
  #   dir.create(path = filepath, showWarnings = TRUE, recursive = TRUE)
  #   filename <- paste0(filepath, "/", c("Class", "Distance"), ".grd")
  #   names(filename) <- c("Class", "Distance")
  # } else if (!canProcessInMemory(r_template, n = length(x@timeseries) + 2)) {
  #   filename <- c(rasterTmpFile("Class"), rasterTmpFile("Distance"))
  # }
  # 
  # if (!is.null(filename)) {
  #   class_b <- writeStart(class_b, filename = filename[1], ...)
  #   distance_b <- writeStart(distance_b, filename = filename[2], ...)
  # } else {
  #   class_vv <- matrix(class_b, ncol = nlayers(class_b))
  #   distance_vv <- matrix(distance_b, ncol = nlayers(distance_b))
  # }
  # 
  # bs <- blockSize(x@timeseries[[1]])
  # bs$array_rows <- cumsum(c(1, bs$nrows * class_vv@ncols))
  # pb <- pbCreate(bs$n, ...)
  # 
  # for(k in 1:bs$n){
  #   # Get raster data
  #   v <- lapply(x@timeseries, getValues, row = bs$row[k], nrows = bs$nrows[k])
  #   
  # }
  
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
    
    levels = c(seq_along(patterns.labels), fill)
    labels = c(patterns.labels, "unclassified")
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
