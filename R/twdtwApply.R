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
#' @param x an object of class twdtw*. This is the target time series. 
#' Usually, it is a set of unclassified time series. 
#'
#' @param y an object of class \link[dtwSat]{twdtwTimeSeries}. 
#' The temporal patterns. 
#' 
#' @param ... arguments to pass to specifique methods for each twdtw* signature 
#' and other arguments for parallel processing to pass to \code{\link[parallel]{mclapply}}.
#'
#' @param resample resample the patterns to have the same length. Default is FALSE.
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
#' Default is FALSE. For \code{\link[dtwSat]{plot-method}} methods use \code{keep=TRUE}.
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
#' @references 
#' [1] Maus  V,  Camara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, vol.PP, no.99, pp.1-11.
#' @references 
#' [2] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' @references 
#' [3] Muller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{twdtw} through the argument \code{weight.fun}. This will 
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
          def = function(x, y, resample=FALSE, length=NULL, weight.fun=NULL, 
                  dist.method="Euclidean", step.matrix = symmetric1, n=NULL, 
                  span=NULL, min.length=0.5, theta = 0.5, keep=FALSE, ...) standardGeneric("twdtwApply"))


#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwTimeSeries 
#' @examples
#' # Applying TWDTW analysis to objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(timeseries=example_ts.list)
#' patterns = twdtwTimeSeries(timeseries=patterns.list, labels=names(patterns.list))
#' matches = twdtwApply(x=ts, y=patterns)
#' matches
#'
#' @export
setMethod(f = "twdtwApply", "twdtwTimeSeries",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep, ...){
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
                  twdtwApply.twdtwTimeSeries(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep, ...)
           })

           
           
#' @rdname twdtwApply 
#' @aliases twdtwApply-twdtwRaster
#' @examples
#' # Applying TWDTW analysis to objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(timeseries=example_ts.list)
#' patterns = twdtwTimeSeries(timeseries=patterns.list, labels=names(patterns.list))
#' matches = twdtwApply(x=ts, y=patterns)
#' matches 
#' 
#' @export
setMethod(f = "twdtwApply", "twdtwRaster",
          def = function(x, y, resample, length, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep, ...){
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
                  twdtwApply.twdtwTimeSeries(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep, ...)
           })
           
           
           
twdtwApply.twdtwTimeSeries = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep, 
mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 1L), mc.cleanup = TRUE){
    if(mc.cores>1){
      res = do.call("joinAlignments", mclapply(as.list(x), mc.preschedule=mc.preschedule, mc.set.seed=mc.set.seed, mc.silent=mc.silent, 
                                       mc.cores=mc.cores, mc.cleanup=mc.cleanup,
                                       FUN = .twdtw, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep))
    } else {
      res = do.call("joinAlignments", lapply(as.list(x), 
                                       FUN = .twdtw, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, keep))
    }
    res
}

