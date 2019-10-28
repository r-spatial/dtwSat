#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwTimeSeries 
#' 
#' @examples 
#' 
#' # Read time series to be labelled from csv
#' ts_file <- system.file("reduce_time/ts_MODIS13Q1.csv", package = "dtwSat")
#' ts <- read.csv(ts_file, stringsAsFactors = FALSE)
#' 
#' # Read labelled temporal patterns from csv
#' pattern_files <- dir(system.file("reduce_time/patterns", package = "dtwSat"), full.names = TRUE)
#' patterns <- lapply(pattern_files, read.csv, stringsAsFactors = FALSE)
#' 
#' # Label time series 
#' ts_class <- .twdtw_reduce_time(x = ts, y = patterns, weight.fun = logisticWeight(-0.1, 50), 
#'                                from = "2009-09-01", to = "2017-09-01", by = "6 month")
#' ts_class
#' 
#' ts_class$pattern_name <- basename(pattern_files)[ts_class$label]
#' ts_class
#' 
#' @export
.twdtw_reduce_time = function(x, 
                              y, 
                              weight.fun = NULL,
                              dist.method = "Euclidean",
                              step.matrix = symmetric1,
                              from = NULL, 
                              to = NULL, 
                              by = NULL, 
                              overlap = .5, 
                              fill = 255){

  # Split time series from dates 
  px <- x[,names(x)!="date",drop=FALSE] 
  tx <- as.Date(x$date) 

  ## Basic use case
  
  # Comput TWDTW alignments for all patterns   
  aligs <- lapply(seq_along(y), function(l){
    
    # Split pattern time series from dates 
    py <- y[[l]][,names(y[[l]])!="date",drop=FALSE] 
    ty <- as.Date(y[[l]]$date)
    
    # Match bands and remove bands that are not in both time series 
    names(py) <- tolower(names(py))
    names(px) <- tolower(names(px))
    px <- px[,names(py),drop=FALSE]
    py <- py[,names(px),drop=FALSE]
    
    # Compute local cost matrix 
    cm <- proxy::dist(py, px, method = dist.method)
    
    if(!is.null(weight.fun)){
      # Get day of the year for pattern and time series 
      doyy <- as.numeric(format(ty, "%j")) 
      doyx <- as.numeric(format(tx, "%j")) 
      
      # Compute time-weght matrix 
      w <- weight.fun(.g(proxy::dist(doyy, doyx, method = dist.method)))
      
      # Apply time-weight to local cost matrix 
      cm <- w * cm
    }
    # Compute accumulated DTW cost matrix 
    internals <- .computecost(cm = cm, step.matrix = step.matrix)
    
    # Find all low cost candidates 
    a <- internals$startingMatrix[internals$N,1:internals$M]
    d <- internals$costMatrix[internals$N,1:internals$M]
    candidates   <- data.frame(a, d)
    candidates   <- candidates[candidates$d==ave(candidates$d, candidates$a, FUN=min),,drop=FALSE]
    candidates$b <- as.numeric(row.names(candidates))
    
    # Order maches by minimum TWDTW distance 
    I <- order(candidates$d)
    if(length(I)<1) return(NULL)

    # Build alignments table table 
    res <- data.frame(
      Alig.N     = I,
      from       = tx[candidates$a[I]],
      to         = tx[candidates$b[I]],
      distance   = candidates$d[I],
      label      = l
    )
    return(res)
    
  })
  
  # Bind rows 
  aligs <- data.table::rbindlist(aligs)

  # Create classification intervals 
  breaks <- seq(as.Date(from), as.Date(to), by = by)
  # Find best macthes for the intervals 
  best_matches <- .bestmatches(
    x = list(aligs[order(aligs$from, aligs$from),,drop=FALSE]), 
    m = length(y), 
    n = length(breaks) - 1, 
    levels = seq_along(y), 
    breaks = breaks, 
    overlap = overlap,
    fill = 99999)$IM

  # Build output 
  out <- as.data.frame(best_matches[,c(1,3),drop=FALSE])
  names(out) <- c("label", "Alig.N")
  out$from <- breaks[-length(breaks)]
  out$to <- breaks[-1]
  out <- merge(out, aligs[, c("label", "Alig.N", "distance"),drop=FALSE], by.x = c("label", "Alig.N"), by.y = c("label", "Alig.N"), all.x = TRUE)
  out <- out[order(out$from), names(out)!="Alig.N"]
  if(any(out$label==0)) out[out$label==0,,drop=FALSE]$label <- fill
  return(out)
}

