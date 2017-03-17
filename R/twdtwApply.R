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
#' a set of time series [1].
#' 
#' @inheritParams twdtwClassify
#' 
#' @param x an object of class twdtw*. This is the target time series. 
#' Usually, it is a set of unclassified time series. 
#'
#' @param y an object of class \link[dtwSat]{twdtwTimeSeries}. 
#' The temporal patterns. 
#' 
#' @param ... arguments to pass to \code{\link[raster]{writeRaster}}
#'
#' @param resample resample the patterns to have the same length. Default is TRUE.
#' See \link[dtwSat]{resampleTimeSeries} for details.
#' 
#' @param length An integer. Patterns length used with \code{patterns.length}. 
#' If not declared the length of the output patterns will be the length of 
#' the longest pattern.
#'  
#' @param weight.fun A function. Any function that receive and performs a 
#' computation on a matrix. The function receives a matrix of time differences 
#' in days and returns a matrix of time-weights. If not declared the time-weight 
#' is zero. In this case the function runs the standard version of the dynamic 
#' time warping. See details. 
#' 
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' see \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' 
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw} [2].
#' 
#' @param n An integer. The maximun number of matches to perform. 
#' NULL will return all matches.
#' 
#' @param theta numeric between 0 and 1. The weight of the time 
#' for the TWDTW computation. Use \code{theta=0} to cancel the time-weight, 
#' \emph{i.e.} to run the original DTW algorithm. Default is 0.5, meaning that 
#' the time has the same weight as the curve shape in the TWDTW analysis.
#' 
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE. For plot methods use \code{keep=TRUE}.
#' 
#' @param span A number. Span between two matches, \emph{i.e.} the minimum  
#' interval between two matches, for details see [3]. If not declared it removes
#' all overlapping matches of the same pattern. To include overlapping matches 
#' of the same pattern use \code{span=0}.
#' 
#' @param min.length A number between 0 an 1. This argument removes the over fittings.
#' Minimum length after warping. Percentage of the original pattern length. Default is 0.5, 
#' meaning that the matching cannot be shorter than half of the pattern length.
#' 
#' @param filepath A character. The path to save the raster with results. If not informed the 
#' function saves in the current work directory. 
#' 
#' @param chunk.size An integer. Set the number of cells for each block, 
#' see \code{\link[raster]{blockSize}} for details.  
#' 
#' @param silent if set to TRUE then hide chunk processing message. Default is FALSE. 
#'
#' @references 
#' [1] Maus  V,  Camara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. IEEE Journal of Selected Topics in Applied Earth Observations and Remote 
#' Sensing, vol.9, no.8, pp.3729-3739.
#' @references 
#' [2] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' @references 
#' [3] Muller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{twdtwApply} through the argument \code{weight.fun}. This will 
#' add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful to analyse time series with phenological cycles
#' of vegetation that are usually bound to seasons. In previous studies by [1] the 
#' logistic weight had better results than the linear for land cover classification. 
#' See [1] for details about the method. 
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
          def = function(x, y, resample=TRUE, length=NULL, weight.fun=NULL, 
                dist.method="Euclidean", step.matrix = symmetric1, n=NULL, 
                span=NULL, min.length=0, theta = 0.5, ...) standardGeneric("twdtwApply"))


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
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep=FALSE, ...){
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
                  twdtwApply.twdtwTimeSeries(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep)
           })

twdtwApply.twdtwTimeSeries = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep){
    res = lapply(as.list(x), FUN = .twdtw, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep)
    new("twdtwMatches", timeseries=x, patterns=y, alignments=res)
}      


           
#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwRaster
#' @examples
#' \dontrun{
#' # Run TWDTW analysis for raster time series 
#' patt = MOD13Q1.MT.yearly.patterns
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
#'                     by="12 month")
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' 
#' r_twdtw = twdtwApply(x=rts, y=patt, weight.fun=log_fun, breaks=time_interval, 
#'           filepath="~/test_twdtw", overwrite=TRUE, format="GTiff")
#'
#' plot(r_twdtw, type="distance")
#' 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' plot(r_lucc)
#' 
#' plot(r_lucc, type="distance")
#' 
#' }
#' @export
setMethod(f = "twdtwApply", "twdtwRaster",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                        breaks=NULL, from=NULL, to=NULL, by=NULL, overlap=0.5, chunk.size=1000, filepath=NULL, silent=FALSE, ...){
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
                  twdtwApply.twdtwRaster(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                                          breaks, overlap, chunk.size, filepath, silent, ...)
           })
           
           
twdtwApply.twdtwRaster = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                                  breaks, overlap, chunk.size, filepath, silent, ...){
     
    # Set blocks Multi-thread parameters
    minblocks = round(nrow(x)*ncol(x) / chunk.size)
    blocks = blockSize(x@timeseries[[1]], minblocks = minblocks)
    threads = seq(1, blocks$n)
    
    # Match raster bands to pattern bands
    raster_bands = coverages(x)
    pattern_names = names(y@timeseries[[1]])
    matching_bands = pattern_names[pattern_names %in% raster_bands]
    if(length(matching_bands)<1)
      stop(paste0("Attributes (bands) of the raster and patterns do not match"))
    if("doy"%in%coverages(x) & !"doy"%in%matching_bands)
      matching_bands = c("doy", matching_bands)
    x = subset(x, layers=matching_bands)
    raster_bands = coverages(x)
    
    # Set raster levels and labels 
    levels = levels(y)
    names(levels) = levels
    
    # Open raster fiels for results 
    if(is.null(filepath)) filepath = paste0(getwd(), "/twdtw_results")
    r_template = brick(x@timeseries[[1]], nl=length(breaks)-1)
    dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
    filename = paste0(filepath,"/twdtw_distance_from_",levels)
    names(filename) = levels
    # b_files = lapply(filename, function(i) writeStart(r_template, filename=i, format="GTiff", overwrite=TRUE))
    b_files = lapply(filename, function(i) writeStart(r_template, filename=i, ...))

    # Get time line 
    timeline = as.Date(index(x))
    
    # Raster info 
    proj_str = projection(x)
    coord    = data.frame(coordinates(x))
    names(coord) = c("longitude","latitude")
    
    fun = function(i){

      if(!silent) print(paste0("Procesing chunk ",i,"/",threads[length(threads)]))
      
      # Get twdtwTimeSeries from raster 
      cells = cellFromRow(x@timeseries[[1]], blocks$row[i]:blocks$nrows[i])
      ts = getTimeSeries(x, y = coord[cells,], proj4string = proj_str)

      # Apply TWDTW analysis  
      twdtw_results = twdtwApply(ts, y=y, weight.fun=weight.fun, dist.method=dist.method, 
                                 step.matrix=step.matrix, n=n, span=span, 
                                 min.length=min.length, theta=theta, keep=FALSE)

      # Get best mathces for each point, period, and pattern 
      m = length(levels)
      h = length(breaks)-1
      A = lapply(as.list(twdtw_results), FUN=.lowestDistances, m=m, n=h, levels=levels, breaks=breaks, overlap=overlap, fill=9999)
      
      # Reshape list to array 
      A = sapply(A, matrix, nrow=h, ncol=m, simplify = 'array')
      
      # Write raster files 
      lapply(seq_along(levels), function(l) writeValues(b_files[[levels[l]]], matrix(t(A[,l,]),ncol=h), blocks$row[i]))
    }
    
    # Apply TWDTW analysis 
    out.list = lapply(threads, FUN = fun)

    # Close raster files 
    b_files = lapply(b_files, writeStop)
    
    # Create brick list for output 
    out = lapply(sapply(b_files, filename), brick)
    
    # Save output timeline 
    timeline = breaks[-1]
    write(as.character(timeline), file = paste(filepath, "timeline", sep="/"))
    
    new("twdtwRaster", timeseries = out, timeline=timeline, layers = names(out))
}

.lowestDistances = function(x, m, n, levels, breaks, overlap, fill){
  .bestmatches(x, m, n, levels, breaks, overlap, fill)$AM
}

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

