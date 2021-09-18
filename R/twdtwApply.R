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
#   R Package dtwSat - 2016-02-22                             #
#                                                             #
###############################################################

#' @include methods.R
#' @title Apply TWDTW analysis 
#' @name twdtwApply
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the matches between the temporal patterns and 
#' a set of time series \insertCite{Maus:2019}{dtwSat}.
#' 
#' @inheritParams twdtwClassify
#' 
#' @param x An object of class twdtw*. This is the target time series. 
#' Usually, it is a set of unclassified time series. 
#'
#' @param y An object of class \link[dtwSat]{twdtwTimeSeries}. 
#' The temporal patterns. 
#' 
#' @param ... Arguments to pass to \code{\link[raster]{writeRaster}} and 
#' \code{\link[raster]{pbCreate}}
#'
#' @param resample Resample the patterns to have the same length. Default is TRUE.
#' See \link[dtwSat]{resampleTimeSeries} for details.
#' 
#' @param length An integer. Length of patterns used with \code{patterns.length}. 
#' If not declared the length of the output patterns will be the length of 
#' the longest pattern.
#'  
#' @param weight.fun A function. Any function that receives two matrices and 
#' performs a computation on them, returning a single matrix with the same 
#' dimensions. The first matrix is the DTW local cost matrix and the 
#' second a matrix of the time differences in days. The function should return a 
#' matrix of DTW local cost weighted by the time differences. If not declared 
#' the time-weight is zero. In this case the function runs the standard version 
#' of the dynamic time warping. See details. 
#' 
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' see \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' 
#' @param step.matrix See \code{\link[dtw]{stepPattern}} in package \pkg{dtw} 
#' \insertCite{Giorgino:2009}{dtwSat}.
#' 
#' @param n An integer. The maximun number of matches to perform. 
#' NULL will return all matches.
#' 
#' @param keep Preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE. For plot methods use \code{keep=TRUE}.
#' 
#' @param span A number. Span between two matches, \emph{i.e.} the minimum  
#' interval between two matches; for details see \insertCite{Muller:2007}{dtwSat}. 
#' If not declared it removes all overlapping matches of the same pattern. To include 
#' overlapping matches of the same pattern use \code{span=0}.
#' 
#' @param min.length A number between 0 an 1. This argument removes overfittings.
#' Minimum length after warping. Percentage of the original pattern length. Default is 0.5, 
#' meaning that the matching cannot be shorter than half of the pattern length.
#' 
#' @param filepath A character. The path at which to save the raster with results. If not provided the 
#' function saves in the current work directory. 
#' 
#' @references 
#'   \insertAllCited{}
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{twdtwApply} through the argument \code{weight.fun}. This will 
#' add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful for analyzing time series with phenological cycles
#' of vegetation that are usually bound to seasons. In previous studies by 
#' \insertCite{Maus:2016}{dtwSat} the logistic weight had better results than the 
#' linear for land cover classification. 
#' See \insertCite{Maus:2016,Maus:2019}{dtwSat} for details about the method. 
#' 
#' @return An object of class twdtw*.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' \code{\link[dtwSat]{getTimeSeries}}, and 
#' \code{\link[dtwSat]{createPatterns}}
#' 
#' @export  
setGeneric(name = "twdtwApply", 
          def = function(x, y, resample=TRUE, length=NULL, weight.fun=function(phi, psi) phi, 
                dist.method="Euclidean", step.matrix = symmetric1, n=NULL, 
                span=NULL, min.length=0, ...) standardGeneric("twdtwApply"))


#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwTimeSeries 
#' @examples
#' # Applying TWDTW analysis to objects of class twdtwTimeSeries
#' log_fun = logisticWeight(-0.1, 100)
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat1 = twdtwApply(x=ts, y=patt, weight.fun=log_fun)
#' mat1
#' 
#' \dontrun{
#' # Parallel processin
#' require(parallel)
#' mat_list = mclapply(as.list(ts), mc.cores=2, FUN=twdtwApply, y=patt, weight.fun=log_fun)
#' mat2 = twdtwMatches(alignments=mat_list)
#' }
#' @export
setMethod(f = "twdtwApply", "twdtwTimeSeries",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, keep=FALSE, ...){
                  if(!is(y, "twdtwTimeSeries"))
                    stop("y is not of class twdtwTimeSeries")
                  if(!is(step.matrix, "stepPattern"))
                    stop("step.matrix is not of class stepPattern")
                  if(is.null(weight.fun))
                    weight.fun = function(psi) 0 
                  if(!is(weight.fun, "function"))
                    stop("weight.fun is not a function")
                  if(resample)
                    y = resampleTimeSeries(object=y, length=length)
                  twdtwApply.twdtwTimeSeries(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, keep)
           })

twdtwApply.twdtwTimeSeries = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, keep){
    res = lapply(as.list(x), FUN = .twdtw, y, weight.fun, dist.method, step.matrix, n, span, min.length, keep)
    new("twdtwMatches", timeseries=x, patterns=y, alignments=res)
}


           
#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwRaster 
#' @example examples/test_twdtw_raster_analysis.R 
#' @export
setMethod(f = "twdtwApply", "twdtwRaster",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, 
                        breaks=NULL, from=NULL, to=NULL, by=NULL, overlap=0.5, filepath="", ...){
                  if(!is(step.matrix, "stepPattern"))
                    stop("step.matrix is not of class stepPattern")
                  if(is.null(weight.fun))
                    weight.fun = function(psi) 0 
                  if(!is(weight.fun, "function"))
                    stop("weight.fun is not a function")
                  if( overlap < 0 & 1 < overlap )
                    stop("overlap out of range, it must be a number between 0 and 1")
                  if(is.null(breaks))
                    if( !is.null(from) &  !is.null(to) ){
                      breaks = seq(as.Date(from), as.Date(to), by=by)    
                    } else {
                      patt_range = lapply(index(y), range)
                      patt_diff = trunc(sapply(patt_range, diff)/30)+1
                      min_range = which.min(patt_diff)
                      by = patt_diff[[min_range]]
                      from = patt_range[[min_range]][1]
                      to = from 
                      month(to) = month(to) + by
                      year(from) = year(range(index(x))[1])
                      year(to) = year(range(index(x))[2])
                      if(to<from) year(to) = year(to) + 1
                      breaks = seq(from, to, paste(by,"month"))
                    }
                  breaks = as.Date(breaks)
                  if(resample)
                    y = resampleTimeSeries(object=y, length=length)
                  twdtwApply.twdtwRaster(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, 
                                          breaks, overlap, filepath, ...)
           })
           

twdtwApply.twdtwRaster = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, 
                                  breaks, overlap, filepath, ...){
  
  
  # Match raster bands to pattern bands
  raster_bands <- coverages(x)
  pattern_bands <- names(y@timeseries[[1]])
  matching_bands <- pattern_bands[pattern_bands %in% raster_bands]
  
  if(length(matching_bands) < 1)
    stop(paste0("Bands from twdtwRaster do not match the bands from patterns"))
  
  if("doy" %in% coverages(x) & !"doy" %in% matching_bands)
    matching_bands <- c("doy", matching_bands)
  
  x <- subset(x, layers = matching_bands)
  raster_bands <- coverages(x)
  
  # Set raster levels and labels 
  levels <- levels(y)
  names(levels) <- levels
  
  # Create output raster objects 
  r_template <- brick(x@timeseries[[1]], nl = length(breaks) - 1, values = FALSE)
  out <- rep(list(r_template), length(levels))
  names(out) <- names(levels)
  
  filepath <- trim(filepath)
  filename <- NULL
  if (filepath != "") {
    dir.create(path = filepath, showWarnings = TRUE, recursive = TRUE)
    filename <- paste0(filepath, "/", names(out), ".grd")
    names(filename) <- names(out)
  } else if (!canProcessInMemory(r_template, n = length(breaks) + length(x@timeseries) )) {
    filename <- sapply(names(out), rasterTmpFile)
  }
  
  if (!is.null(filename)) {
    out <- lapply(names(out), function(i) writeStart(out[[i]], filename = filename[i], ...))
  } else {
    vv <- lapply(names(out), function(i) matrix(out[[i]], ncol = nlayers(out[[i]])))
    names(vv) <- names(out)
  }
  
  bs <- blockSize(x@timeseries[[1]])
  bs$array_rows <- cumsum(c(1, bs$nrows*out[[1]]@ncols))
  pb <- pbCreate(bs$n, ...)
  
  for(k in 1:bs$n){
    
    # Get raster data
    v <- lapply(x@timeseries, getValues, row = bs$row[k], nrows = bs$nrows[k])
    
    # Create twdtwTimeSeries
    ts <- dtwSat::twdtwTimeSeries(lapply(1:nrow(v[[1]]), function(i){
      # Get time series dates
      if(any(names(v) %in% "doy")){
        timeline <- dtwSat::getDatesFromDOY(year = format(dtwSat::index(x), "%Y"), doy = v[["doy"]][i,])
      } else {
        timeline <- dtwSat::index(x)
      }
      # Tag duplicate dates for removal
      s = !duplicated(timeline)
      # Create multi band time series
      zoo::zoo(sapply(v[!(names(v) %in% "doy")], function(v) v[i, s, drop = FALSE]), timeline[s])
    }))
    
    # Apply TWDTW analysis
    twdtw_results <- dtwSat::twdtwApply(x = ts, y = y, weight.fun = weight.fun, dist.method = dist.method,
                                        step.matrix = step.matrix, n = n, span = span,
                                        min.length = min.length, keep = FALSE)
    
    # Get best matches for each point, period, and pattern
    m <- length(levels)
    h <- length(breaks) - 1
    
    A <- lapply(as.list(twdtw_results), FUN = .lowestDistances, m = m, n = h, 
                levels = levels, breaks = breaks, overlap = overlap, fill = 9999)
    
    B <- as.list(data.frame(do.call("rbind", A)))
    
    names(B) <- levels
    
    B <- lapply(B, function(m) matrix(as.numeric(m), nrow = length(A), ncol = nrow(A[[1]]), byrow = TRUE))
    
    rows <- seq(from = bs$array_rows[k], by = 1, length.out = bs$nrows[k]*out[[1]]@ncols)
    for(l in seq_along(levels)){
      if (!is.null(filename)) {
        writeValues(out[[l]], B[[l]], bs$row[k])
      } else {
        vv[[l]][rows,] <- B[[l]]
      }
    }
    
    pbStep(pb, k)
    
  }
  
  if (!is.null(filename)) {
    out <- lapply(out, writeStop)
  } else {
    out <- lapply(seq_along(levels), function(i) setValues(out[[i]], values = vv[[i]]))
  }
  
  pbClose(pb)
  
  names(out) <- levels
  
  new("twdtwRaster", timeseries = out, timeline = breaks[-1], layers = levels)
  
}

.lowestDistances = function(x, m, n, levels, breaks, overlap, fill){
  .bestmatches(x, m, n, levels, breaks, overlap, fill)$AM
}

# .TlowestDistances = function(x, m, n, levels, breaks, overlap, fill){
#   t(.bestmatches(x, m, n, levels, breaks, overlap, fill)$AM)
# }

# Crop raster time series. Returns a 3D array 
.cropTimeSeries = function(x, r1, r2){
  if(is(x, "RasterBrick") | is(x, "RasterStack")){
    y = extent(x, r1, r2)
    x = crop(x, y)
  } else {
    x = x[r1:r2,,]
  }
  alply(x, c(1,2), as.numeric)
}

# Build zoo time series  
.bulidZoo = function(p, x, timeline){
  # Get time series for each band 
  datasets = lapply(x, function(x) x[p,])
  if(any("doy"==names(datasets))){
    datasets$doy = getDatesFromDOY(doy=datasets$doy, year=format(timeline, "%Y"))
  } else {
    datasets$doy = timeline
  }
    
  idoy = which(names(datasets) %in% c("doy"))
  
  # Remove invalid values 
  k = unlist(lapply(datasets[-idoy], function(x){
    which(x<0|is.na(x))
  }))
  k = c(k, which(duplicated(datasets$doy)))
  if(length(k)>0) datasets = lapply(datasets, function(x) x[-k] )
  
  # Build multi-band zoo object 
  zoo(data.frame(datasets[-idoy]), order.by = datasets$doy)
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

