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
           def = function(x, y, resample = TRUE, length = NULL, weight.fun = NULL, 
                          dist.method = "Euclidean", step.matrix = symmetric1, n = NULL, 
                          span = NULL, min.length = 0, theta = 0.5, ...) standardGeneric("twdtwApplyParallel"))



#' @rdname twdtwApplyParallel 
#' @aliases twdtwApplyParallel-twdtwRaster
# # \dontrun{
#' @examples
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
#' plot(r_lucc)
#' 
#' # Assess classification 
#' twdtw_assess = twdtwAssess(object = r_lucc, y = validation_samples, 
#'                            proj4string = proj_str, conf.int = .95, rm.nosample = TRUE) 
#' twdtw_assess
#' 
#' # Plot assessment 
#' plot(twdtw_assess, type="accuracy")
#' plot(twdtw_assess, type="area")
#' 
#' @export
setMethod(f = "twdtwApplyParallel", "twdtwRaster",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
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
            twdtwApplyParallel.twdtwRaster(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                                   breaks, overlap, filepath, ...)
          })


twdtwApplyParallel.twdtwRaster = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
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
  } else if (!canProcessInMemory(r_template, n = length(breaks) + length(x@timeseries))) {
    filename <- sapply(names(out), rasterTmpFile)
  }
  # filename <- sapply(names(out), rasterTmpFile)
  
  if (!is.null(filename)) {
    out <- lapply(names(out), function(i) writeStart(out[[i]], filename = filename[i], ...))
  } else {
    vv <- lapply(names(out), function(i) matrix(out[[i]], ncol = nlayers(out[[i]])))
    names(vv) <- names(out)
  }
  
  cl <- getCluster()
  on.exit( returnCluster() )
  
  nodes <- length(cl)
  
  bs <- blockSize(x@timeseries[[1]], minblocks = nodes*4)
  bs$array_rows <- cumsum(c(1, bs$nrows*out[[1]]@ncols))
  pb <- pbCreate(bs$n, ...)

  clusterExport(cl = cl, list = c("y", "weight.fun", "dist.method", "step.matrix", "n", "span", 
                                  "min.length", "theta", "breaks", "overlap"), envir = environment())
  
  fun <- function(k){
    
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
                                       min.length = min.length, theta = theta, keep = FALSE)

    # Get best mathces for each point, period, and pattern
    levels <- dtwSat::levels(y)
    m <- length(levels)
    h <- length(breaks)-1

    A <- lapply(as.list(twdtw_results), FUN = .lowestDistances, m = m, n = h, 
                levels = levels, breaks = breaks, overlap = overlap, fill = 9999)
    B <- as.list(data.frame(do.call("rbind", A)))
    names(B) <- levels
    lapply(B, function(m) matrix(as.numeric(m), nrow = length(A), ncol = nrow(A[[1]]), byrow = TRUE))
  }
  
  # get all nodes going
  for (k in 1:nodes) {
    sendCall(cl[[k]], fun, list(k), tag = k)
  }
  
  for (k in 1:bs$n) {
    # receive results from a node
    d <- recvOneData(cl)
    
    # error?
    if (!d$value$success) { 
      stop('cluster error')
    }
    
    # which block is this?
    b <- d$value$tag
    
    if (!is.null(filename)) {
      out <- lapply(seq_along(levels), function(l) writeValues(out[[l]], d$value$value[[l]], bs$row[b]))
    } else {
      rows <- seq(from = bs$array_rows[b], by = 1, length.out = bs$nrows[b]*out[[1]]@ncols)
      for(l in seq_along(levels)){
        vv[[l]][rows,] <- d$value$value[[l]]
      }
    }
    
    # need to send more data?
    ni <- nodes + k
    if (ni <= bs$n) {
      sendCall(cl[[d$node]], fun, list(ni), tag = ni)
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




