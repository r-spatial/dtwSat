#' @include methods.R
#' @title Apply TWDTW analysis to twdtwRaster using parallel processing 
#' @name twdtwApplyParallel 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the matches between the temporal patterns and 
#' a set of time series [1].
#' 
#' @inheritParams twdtwApply
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
#' @return An object of class twdtwRaster.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}}, and 
#' \code{\link[dtwSat]{createPatterns}}
#' 
#' @export  
setGeneric(name = "twdtwApplyParallel", 
           def = function(x, y, resample=TRUE, length=NULL, weight.fun=NULL, 
                          dist.method="Euclidean", step.matrix = symmetric1, n=NULL, 
                          span=NULL, min.length=0, theta = 0.5, ...) standardGeneric("twdtwApplyParallel"))



#' @rdname twdtwApplyParallel 
#' @aliases twdtwApplyParallel-twdtwRaster
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
setMethod(f = "twdtwApplyParallel", "twdtwRaster",
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
            twdtwApplyParallel.twdtwRaster(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                                   breaks, overlap, chunk.size, filepath, silent, ...)
          })


twdtwApplyParallel.twdtwRaster = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
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




# Prallel raster processing 
apply_raster_parallel <- function(x, fun, A, filepath="", ...) {
  
  # Set blocks Multi-thread parameters
  # minblocks <- round(nrow(x) * ncol(x) / chunk.size)
  # blocks <- blockSize(x@timeseries[[1]], minblocks = minblocks)
  # threads <-  seq(1, blocks$n)
  
  # Match raster bands to pattern bands
  raster_bands <- coverages(x)
  pattern_names <- names(y@timeseries[[1]])
  matching_bands <- pattern_names[pattern_names %in% raster_bands]
  
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
  # if(is.null(filepath)) 
  #   filepath = paste0(getwd(), "/twdtw_results")
  r_template <- brick(x@timeseries[[1]], nl = length(breaks) - 1, values = FALSE)
  out <- rep(list(r_template), length(levels))
  names(out) <- names(levels)
  # dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
  # filename = paste0(filepath,"/twdtw_distance_from_",levels)
  # names(filename) = levels
  # b_files = lapply(filename, function(i) writeStart(r_template, filename=i, format="GTiff", overwrite=TRUE))
  # b_files = lapply(filename, function(i) writeStart(r_template, filename=i, ...))
  
  # out <- brick(x@timeseries[[1]], nl=length(breaks)-1, values=FALSE)
  # names(out) <- names(x)
  
  cl <- getCluster()
  on.exit( returnCluster() )
  
  nodes <- length(cl)
  
  bs <- blockSize(x, minblocks = nodes*4)
  bs$array_rows <- cumsum(c(1, bs$nrows*out@ncols))
  pb <- pbCreate(bs$n)
  
  clFun <- function(k){
    
    # Get raster data 
    # v <- lapply(x@timeseries, getValues, row = 1, nrows = 2)
    v <- lapply(x, getValues, row = bs$row[k], nrows = bs$nrows[k])
    
    # Create twdtwTimeSeries 
    ts <- twdtwTimeSeries(lapply(1:nrow(v[[1]]), function(i){
      # Get time series dates 
      if(any(names(v) %in% "doy")){
        timeline <- getDatesFromDOY(year = format(index(x), "%Y"), doy = v[["doy"]][i,])
      } else {
        timeline <- index(x)
      }
      # Tag duplicate dates for removal 
      k = !duplicated(timeline)
      # Create multi band time series 
      zoo(sapply(v[!(names(v) %in% "doy")], function(v) v[i, k, drop = FALSE]), timeline[k])
    }))
    
    # Apply TWDTW analysis  
    twdtw_results = twdtwApply(ts, y=y, weight.fun=weight.fun, dist.method=dist.method, 
                               step.matrix=step.matrix, n=n, span=span, 
                               min.length=min.length, theta=theta, keep=FALSE)
    
    # Get best mathces for each point, period, and pattern 
    m = length(levels)
    h = length(breaks)-1
  
    A <- lapply(as.list(twdtw_results), FUN=.lowestDistances, m=m, n=h, levels=levels, breaks=breaks, overlap=overlap, fill=9999)
    B <- as.list(data.frame(do.call("rbind", A)))
    lapply(B, function(m) matrix(as.numeric(m), nrow = length(A), ncol = nrow(A[[1]]), byrow = TRUE))
  }
  
  # get all nodes going
  for (k in 1:nodes) {
    sendCall(cl[[k]], clFun, list(k, A), tag = k)
  }
  
  filepath <- trim(filepath)
  filename <- NULL
  if (filepath != "") {
    filename <- paste0(filepath, "/", names(out), ".grd")
  } else if (!canProcessInMemory(r_template, n = length(breaks))) {
    filename <- sapply(names(out), rasterTmpFile)
  } 
  
  if (!is.null(filename)) {
    # out <- writeStart(out, filename = filename, ... )
    out <- lapply(names(out), function(i) writeStart(out[[i]], filename = filename[i], ...))
  } else {
    # vv <- matrix(ncol=nrow(out), nrow=ncol(out))
    # vv <- matrix(out, ncol=nlayers(out))
    vv <- lapply(names(out), function(i) matrix(out[[i]], ncol = nlayers(out[[i]])))
  }
  
  for (k in 1:bs$n) {
    # receive results from a node
    d <- recvOneData(cl)
    
    # error?
    if (! d$value$success) { 
      stop('cluster error')
    }
    
    # which block is this?
    b <- d$value$tag
    cat('received block: ',b," / ",bs$n,'\n'); flush.console();
    
    if (!is.null(filename)) {
      # out <- writeValues(out, d$value$value, bs$row[b])
      out <- lapply(seq_along(levels), function(l) writeValues(out[[l]], d$value$value[[l]], bs$row[b]))
    } else {
      # cols <- bs$row[b]:(bs$row[b] + bs$nrows[b]-1)	
      # vv[,cols] <- matrix(d$value$value, nrow=out@ncols)
      rows <- seq(from = bs$array_rows[b], by = 1, length.out = bs$nrows[b]*out@ncols)
      vv <- lapply(seq_along(levels), function(l) vv[[l]][rows,] <- d$value$value[[l]])
    }
    
    # need to send more data?
    ni <- nodes + k
    if (ni <= bs$n) {
      sendCall(cl[[d$node]], clFun, list(ni, A), tag = ni)
    }	
    pbStep(pb)
  }
  if (!is.null(filename)) {
    # out <- writeStop(out)
    out <- lapply(out, writeStop)
  } else {
    # out <- setValues(out, as.vector(vv))
    out <- lapply(names(out), function(i) setValues(out[[i]], values = as.vector(vv[[i]])))
  }
  pbClose(pb)
  
  return(out)
}

.lowestDistances2 = function(x, m, n, levels, breaks, overlap, fill){
  t(.bestmatches(x, m, n, levels, breaks, overlap, fill)$AM)
}
