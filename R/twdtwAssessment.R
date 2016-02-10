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
#   R Package dtwSat - 2016-02-10                             #
#                                                             #
###############################################################


#' @title Assess TWDTW analyis using validation samples  
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function classifies the set of samples time series 
#' 
#' @param x A \code{\link[base]{list}} of three attributes: the trained \code{pattrens},
#' the validation dataset \code{x}, and the reference for the validation dataset \code{ref}.
#' See \code{\link[dtwSat]{splitDataset}} for details. 
#' 
#' @param overlap A number between 0 and 1. The minimum overlapping between one match and the 
#' interval of classification. Default is 0.3, i.e. an overlap minimum of 30\%. 
#' See \code{\link[dtwSat]{classifyIntervals}} for details. 
#' 
#' @param threshold A number. The TWDTW dissimilarity threshold, \emph{i.e.} the maximum TWDTW cost 
#' for consideration in the classification. Default is Inf.
#' See \code{\link[dtwSat]{classifyIntervals}} for details. 
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{twdtw}}.
#' 
#' 
#' @docType methods
#' @return A \code{\link[base]{data.frame}} of two attributes \code{Reference} 
#' and \code{Predicted}.
#' 
#' @seealso 
#' \code{\link[dtwSat]{splitDataset}}
#' 
#' @examples
#' 
#' # set.seed(1)
#' # field_samples = read.csv(system.file("lucc_MT/samples.csv", package="dtwSat")) 
#' # proj_str = scan(file = system.file("lucc_MT/samples_projection", package="dtwSat"), 
#' #                  what = "character")
#' # reference = as.character(field_samples[["class"]])
#' # load(system.file("lucc_MT/ts_list.RData", package="dtwSat"))
#' 
#' # x = splitted_dataset = splitDataset(timeseries = ts_list, ref = reference, 
#' #        times=2, p=0.1, mc.cores=1, freq=8, from="2007-09-01", to="2008-09-01", 
#' #        formula = y ~ s(time, bs="cc"))
#' 
#' # assess_table = lapply(x, FUN=twdtwAssessment, overlap = 0.5, 
#' #             weight.fun = logisticWeight(alpha=-0.1, beta=50))
#' 
#' @export
twdtwAssessment = function(x, mc.cores=1, overlap=0.5, threshold=Inf, ...){
  # Apply TWDTW for each time series 
  twdtw_results = mclapply(x$x, mc.cores = mc.cores, FUN = twdtw, patterns = x$patterns, ...)
  # Classify time interval for each pixel 
  pred = sapply(twdtw_results, function(s){
    ts = getTimeSeries(s)
    classifyIntervals(x = s, from=min(index(ts)), to=max(index(ts)), by=max(index(ts))-min(index(ts)), overlap=overlap, threshold=threshold)$pattern
  })
  res = data.frame(Reference=x$ref, Predicted=pred)
  res
}

