#' @include methods.R
#' @title Minimalist version of TWDTW apply
#' @name twdtwReduceTime
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' @rdname twdtwReduceTime 
#' 
#' @description This function is a minimalist implementation of 
#' \link[dtwSat]{twdtwApply} that is in average 3x faster. The time weight function 
#' is coded in Fortran. It does not keep any intermediate data. 
#' It performs a multidimensional TWDTW analysis 
#' \insertCite{Maus:2019}{dtwSat} and retrieves only the best matches between 
#' the unclassified time series and the patterns for each defined time interval.
#' 
#' @inheritParams twdtwApply
#' @inheritParams twdtwClassify
#' 
#' @param x a data.frame with the target time series. Usually, it is an 
#' unclassified time series. It must contain two or more columns, one column 
#' called \code{date} with dates in the format "YYYY-MM-DD". The other columns 
#' can have any names (e.g., red, blue, nir, evi, ndvi) as long as they match 
#' the column names in the temporal patterns \code{y}. 
#' 
#' @param y a list of data.frame objects similar to \code{x}. 
#' The temporal patterns used to classify the time series in \code{x}. 
#' 
#' @param fill An integer to fill the classification gaps. Default 255.
#' 
#' @examples 
#' \dontrun{
#' 
#' library(dtwSat)
#' from = "2009-09-01"
#' to = "2017-09-01"
#' by = "12 month"
#'
#' # S4 objects for original implementation 
#' tw_patt = readRDS(system.file("lucc_MT/patterns/patt.rds", package = "dtwSat"))
#' tw_ts = twdtwTimeSeries(MOD13Q1.ts) 
#' 
#' # Table from csv for minimalist version 
#' mn_patt <- lapply(dir(system.file("lucc_MT/patterns", package = "dtwSat"), 
#'   pattern = ".csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE)
#' mn_ts <- read.csv(system.file("reduce_time/ts_MODIS13Q1.csv", package = "dtwSat"), 
#'   stringsAsFactors = FALSE)
#' 
#' # Benchtmark 
#' rbenchmark::benchmark(
#'   original = twdtwClassify(twdtwApply(x = tw_ts, y = tw_patt, weight.fun = log_fun), 
#'                                       from = from, to = to, by = by)[[1]],
#'   minimalist = twdtwReduceTime(x = mn_ts, y = mn_patt, 
#'                                       from = from, to = to, by = by)  
#'  )
#' }
#' 
#' @export
twdtwReduceTime = function(x, 
                           y, 
                           alpha = -0.1,
                           beta = 50,
                           dist.method = "Euclidean",
                           step.matrix = symmetric1,
                           from = NULL, 
                           to = NULL, 
                           by = NULL, 
                           breaks = NULL,
                           overlap = .5, 
                           fill = 255){

  # Split time series from dates 
  px <- x[,names(x)!="date",drop=FALSE] 
  tx <- as.Date(x$date) 

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
    #cm <- proxy::dist(py, px, method = dist.method)
    
    # Get day of the year for pattern and time series 
    doyy <- as.numeric(format(ty, "%j")) 
    doyx <- as.numeric(format(tx, "%j")) 
    
    # Compute accumulated DTW cost matrix 
    xm = na.omit(cbind(doyx, as.matrix(px)))
    ym = na.omit(cbind(doyy, as.matrix(py)))
    internals = .fast_twdtw(xm, ym, alpha, beta, step.matrix)
    
    # Find all low cost candidates 
    a <- internals$startingMatrix[internals$N-1,1:internals$M]
    d <- internals$costMatrix[internals$N-1,1:internals$M]
    # rbenchmark::benchmark(
    #   df = data.frame(a, d),
    #   mt = matrix(c(a, d, 1:internals$M), ncol = 3, byrow = F), replications = 100)
    # candidates   <- data.frame(a, d)
    candidates   <- matrix(c(a, d, 1:internals$M, 1:internals$M, rep(l, internals$M)), ncol = 5, byrow = F)
    # candidates   <- candidates[candidates$d==ave(candidates$d, candidates$a, FUN=min),,drop=FALSE]
    candidates   <- candidates[candidates[,2]==ave(candidates[,2], candidates[,1], FUN=min),,drop=FALSE]
    # candidates$b <- as.numeric(row.names(candidates))
    # candidates[,2] <- as.numeric(row.names(candidates))
    candidates <- candidates[!is.na(candidates[,1]),,drop=FALSE]
    
    # Order maches by minimum TWDTW distance 
    # I <- order(candidates$d)
    I <- order(candidates[,3])
    if(length(I)<1) return(NULL)

    # Build alignments table table 
    candidates[,4] <- I
    # res <- data.frame(
    #   Alig.N     = I,
    #   from       = tx[candidates[I,1]],
    #   to         = tx[candidates[I,3]],
    #   distance   = candidates[I,2],
    #   label      = l
    # )
    return(candidates)
    
  })
  
  # Bind rows 
  aligs <- do.call("rbind", aligs)
  
  # Create classification intervals 
  if(is.null(breaks)){
    breaks <- seq(as.Date(from), as.Date(to), by = by) 
  }
  
  # Find best macthes for the intervals 
  best_matches <- .bestmatches2(
    x = aligs[order(aligs[,1], aligs[,2]),,drop=FALSE], 
    tx = tx,
    m = length(y), 
    n = length(breaks) - 1, 
    levels = seq_along(y), 
    breaks = breaks, 
    overlap = overlap,
    fill = 99999)[c("IM", "DB")]

  # Build output 
  out <- as.data.frame(best_matches$IM[,1,drop=FALSE])
  names(out) <- c("label")
  out$from <- breaks[-length(breaks)]
  out$to <- breaks[-1]
  out$distance <- best_matches$DB
  # out <- merge(out, aligs[, c("label", "Alig.N", "distance"),drop=FALSE], by.x = c("label", "Alig.N"), by.y = c("label", "Alig.N"), all.x = TRUE)
  # out <- out[order(out$from), names(out)!="Alig.N"]
  if(any(out$label==0)) out[out$label==0,]$label <- fill
  return(out)
}

# @useDynLib dtwSat computecost
.fast_twdtw = function(xm, ym, alpha, beta, step.matrix){
  
  #  cm = rbind(0, cm)
  n = nrow(ym)
  m = nrow(xm)
  d = ncol(ym)
  
  if(is.loaded("twdtw", PACKAGE = "dtwSat", type = "Fortran")){
    out = .Fortran(twdtw, 
                   XM = matrix(as.double(xm), m, d),
                   YM = matrix(as.double(ym), n, d),
                   CM = matrix(as.double(0), n+1, m),
                   DM = matrix(as.integer(0), n+1, m),
                   VM = matrix(as.integer(0), n+1, m),
                   SM = matrix(as.integer(step.matrix), nrow(step.matrix), ncol(step.matrix)),
                   N  = as.integer(n),
                   M  = as.integer(m),
                   D  = as.integer(d),
                   NS = as.integer(nrow(step.matrix)),
                   TW = as.double(c(alpha, beta))
                   # WB = as.double(beta)
                   )
  } else {
    stop("Fortran twdtw lib is not loaded")
  }
  # sqrt(sum((ym[1,-1] - xm[1,-1])*(ym[1,-1] - xm[1,-1])))
  res = list()
  res$costMatrix = out$CM[-1,]
  res$directionMatrix = out$DM[-1,]
  res$startingMatrix = out$VM[-1,]
  res$stepPattern = step.matrix
  res$N = n
  res$M = m
  res
}


# @useDynLib dtwSat bestmatches
.bestmatches2 = function(x, tx, m, n, levels, breaks, overlap, fill=9999){
  if(is.loaded("bestmatches", PACKAGE = "dtwSat", type = "Fortran")){
    if(length(x[,2])<1){
      res = list(
        XM = matrix(as.integer(c(as.numeric(tx[x[,1]]), as.numeric(tx[x[,3]]))), ncol = 2),
        AM = matrix(as.double(fill), nrow = n, ncol = m), 
        DM = as.double(x[,2]),
        DP = as.integer(as.numeric(breaks)),
        X  = as.integer(match(x[,5], levels)),
        IM = matrix(as.integer(0), nrow = n, ncol = 3),
        DB = as.double(x[,2]),
        A  = as.integer(x[,4]),
        K  = as.integer(length(x[,4])),
        P  = as.integer(length(breaks)),
        L  = as.integer(length(levels)),
        OV = as.double(overlap))
    } else {
      res = try(.Fortran(bestmatches, 
                         XM = matrix(as.integer(c(as.numeric(tx[x[,1]]), as.numeric(tx[x[,3]]))), ncol = 2),
                         AM = matrix(as.double(fill), nrow = n, ncol = m), 
                         DM = as.double(x[,2]),
                         DP  = as.integer(as.numeric(breaks)),
                         X  = as.integer(match(x[,5], levels)),
                         IM = matrix(as.integer(0), nrow = n, ncol = 3),
                         DB = as.double(rep(0, n)),
                         A  = as.integer(x[,4]),
                         K  = as.integer(length(x[,4])),
                         P  = as.integer(length(breaks)),
                         L  = as.integer(length(levels)),
                         OV = as.double(overlap))) 
    }
  } else {
    stop("Fortran bestmatches lib is not loaded")
  }
  if(is(res, "try-error")){
    res = list(
      XM = matrix(as.integer(c(as.numeric(tx[x[,1]]), as.numeric(tx[x[,3]]))), ncol = 2),
      AM = matrix(as.double(fill), nrow = n, ncol = m), 
      DM = as.double(x[,2]),
      DP = as.integer(as.numeric(breaks)),
      X  = as.integer(match(x[,5], levels)),
      IM = matrix(as.integer(0), nrow = n, ncol = 3),
      DB = as.double(x[,2]),
      A  = as.integer(x[,4]),
      K  = as.integer(length(x[,4])),
      P  = as.integer(length(breaks)),
      L  = as.integer(length(levels)),
      OV = as.double(overlap)
    )
  } 
  res
}


