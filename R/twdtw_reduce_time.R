#' @include methods.R
#' @title Faster version of TWDTW apply
#' @name twdtwReduceTime
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' @rdname twdtwReduceTime 
#' 
#' @description This function is a faster implementation of 
#' \link[dtwSat]{twdtwApply} that is in average 4x faster. The time weight function 
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
#' @param time.window logical. TRUE will constrain the TWDTW computation to the 
#' value of the parameter \code{beta} defined in the logistic weight function. 
#' Default is FALSE. 
#' 
#' @param fill An integer to fill the classification gaps.
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
#' # Table from csv for faster version 
#' mn_patt <- lapply(dir(system.file("lucc_MT/patterns", package = "dtwSat"), 
#'   pattern = ".csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE)
#' mn_ts <- read.csv(system.file("reduce_time/ts_MODIS13Q1.csv", package = "dtwSat"), 
#'   stringsAsFactors = FALSE)
#' 
#' # Benchtmark 
#' rbenchmark::benchmark(
#'   legacy_twdtw = twdtwClassify(twdtwApply(x = tw_ts, y = tw_patt, weight.fun = log_fun), 
#'                                       from = from, to = to, by = by)[[1]],
#'   fast_twdtw = twdtwReduceTime(x = mn_ts, y = mn_patt, rom = from, to = to, by = by)  
#'  )
#' }
#' 
#' @export
twdtwReduceTime = function(x,
                           y,
                           alpha = -0.1,
                           beta = 50,
                           time.window = FALSE,
                           dist.method = "Euclidean",
                           step.matrix = symmetric1,
                           from = NULL,
                           to = NULL,
                           by = NULL,
                           breaks = NULL,
                           overlap = .5,
                           fill = length(y) + 1,
                           keep = FALSE, ...){

  # Split time series from dates
  px <- x[,names(x)!="date",drop=FALSE]
  tx <- as.Date(x$date)

  # Compute TWDTW alignments for all patterns   
  twdtw_data <- lapply(seq_along(y), function(l){

    # Split pattern time series from dates 
    py <- y[[l]][,names(y[[l]])!="date",drop=FALSE] 
    ty <- as.Date(y[[l]]$date)
    
    # Match bands and remove bands that are not in both time series 
    names(py) <- tolower(names(py))
    names(px) <- tolower(names(px))
    px <- px[,names(py),drop=FALSE]
    py <- py[,names(px),drop=FALSE]
    
    # Get day of the year for pattern and time series 
    doyy <- as.numeric(format(ty, "%j")) 
    doyx <- as.numeric(format(tx, "%j")) 
    
    # Compute accumulated DTW cost matrix 
    xm = na.omit(cbind(doyx, as.matrix(px)))
    ym = na.omit(cbind(doyy, as.matrix(py)))
    internals = .fast_twdtw(xm, ym, alpha, beta, step.matrix, time.window)
    
    # Find all low cost candidates 
    b <- internals$JB[internals$JB!=0]
    a <- internals$VM[-1,][internals$N,b]
    d <- internals$CM[-1,][internals$N,b]
    # View(internals$VM[-1,])
    # CM <- internals$CM[-1,]; CM[CM>10000] <- NA; CM |> t() |> image()
    candidates <- matrix(c(a, d, b, b, rep(l, length(b))), ncol = 5, byrow = F)
    
    # Order matches by minimum TWDTW distance 
    I <- order(candidates[,3])
    if(length(I)<1) return(NULL)

    # Build alignments table table 
    candidates[,4] <- I
    
    if(!keep){
      internals = NULL
    } else {
      internals$tsDates <- tx
      internals$patternDates <- ty
    }
    
    return(list(candidates = candidates, internals = internals))
    
  })
  
  aligs <- do.call("rbind", lapply(twdtw_data, function(x) x$candidates))
  il <- order(aligs[,1], aligs[,2])
  
  # Create classification intervals 
  # if(is.null(breaks)){
  #   breaks <- seq(as.Date(from), as.Date(to), by = by) 
  # }
  if(is.null(breaks))
    if( !is.null(from) & !is.null(to) ){
      breaks = seq(as.Date(from), as.Date(to), by = by)
    } else {
      patt_range = lapply(y, function(yy) range(yy$date))
      patt_diff = trunc(sapply(patt_range, diff)/30)+1
      min_range = which.min(patt_diff)
      by = patt_diff[[min_range]]
      from = patt_range[[min_range]][1]
      to = from 
      month(to) = month(to) + by
      year(from) = year(range(x$date)[1])
      year(to) = year(range(x$date)[2])
      if(to<from) year(to) = year(to) + 1
      breaks = seq(from, to, paste(by,"month"))
      breaks = as.Date(breaks)
    }
  
  # Find best macthes for the intervals 
  best_matches <- .bestmatches2(
    x = aligs[il,,drop=FALSE], 
    tx = tx,
    m = length(y), 
    n = length(breaks) - 1, 
    levels = seq_along(y), 
    breaks = breaks, 
    overlap = overlap,
    fill = fill)

  # Build output 
  out <- list(
    label = best_matches$IM[,1],
    from = breaks[-length(breaks)],
    to = breaks[-1],
    distance = best_matches$DB)
  
  # Bind rows 
  if(keep){
    aligs <- cbind(aligs, rep(0, nrow(aligs)))
    aligs[il[best_matches$IM[best_matches$IM[,2]!=fill,2]], 6] <- 1
    out$internals <- list(alignments = aligs, internals = lapply(twdtw_data, function(x) x$internals))
  }
  
  return(out)
}

# @useDynLib dtwSat computecost
.fast_twdtw = function(xm, ym, alpha, beta, step.matrix, wc){
  
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
                   TW = as.double(c(alpha, beta)),
                   LB = as.logical(wc),
                   JB = as.integer(rep(0, n))
                   )
  } else {
    stop("Fortran twdtw lib is not loaded")
  }
  out
}


# @useDynLib dtwSat bestmatches
.bestmatches2 = function(x, tx, m, n, levels, breaks, overlap, fill){
  if(is.loaded("bestmatches", PACKAGE = "dtwSat", type = "Fortran")){
    if(length(x[,2])<1){
      res = list(
        XM = matrix(as.integer(c(as.numeric(tx[x[,1]]), as.numeric(tx[x[,3]]))), ncol = 2),
        AM = matrix(as.double(.Machine$double.xmax), nrow = n, ncol = m), 
        DM = as.double(x[,2]),
        DP = as.integer(as.numeric(breaks)),
        X  = as.integer(match(x[,5], levels)),
        IM = matrix(as.integer(fill), nrow = n, ncol = 3),
        DB = as.double(x[,2]),
        A  = as.integer(x[,4]),
        K  = as.integer(length(x[,4])),
        P  = as.integer(length(breaks)),
        L  = as.integer(length(levels)),
        OV = as.double(overlap))
    } else {
      res = try(.Fortran(bestmatches, 
                         XM = matrix(as.integer(c(as.numeric(tx[x[,1]]), as.numeric(tx[x[,3]]))), ncol = 2),
                         AM = matrix(as.double(.Machine$double.xmax), nrow = n, ncol = m), 
                         DM = as.double(x[,2]),
                         DP  = as.integer(as.numeric(breaks)),
                         X  = as.integer(match(x[,5], levels)),
                         IM = matrix(as.integer(fill), nrow = n, ncol = 3),
                         DB = as.double(rep(.Machine$double.xmax, n)),
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
      AM = matrix(as.double(.Machine$double.xmax), nrow = n, ncol = m), 
      DM = as.double(x[,2]),
      DP = as.integer(as.numeric(breaks)),
      X  = as.integer(match(x[,5], levels)),
      IM = matrix(as.integer(fill), nrow = n, ncol = 3),
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

