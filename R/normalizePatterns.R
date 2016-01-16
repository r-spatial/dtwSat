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
#   R Package dtwSat - 2016-16-01                             #
#                                                             #
###############################################################


#' @title Normalize patterns length 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function normalizes the length of a pattern or 
#' list of patterns
#' 
#' @param ... \link[zoo]{zoo} objects
#' @param patterns a list of \link[zoo]{zoo} objects
#' @param patterns.length An integer. Default is the length of the longest pattern
#' @docType methods
#' @return A \link[base]{list} of \link[zoo]{zoo} objects with normalized 
#' length 
#'
#' @examples
#' new.patterns.list = normalizePatterns(patterns = patterns.list, patterns.length = 23)
#' sapply(patterns.list, nrow)
#' sapply(new.patterns.list, nrow)
#' 
#' @export
#' 
normalizePatterns = function(..., patterns = list(...), patterns.length=NULL){
  if(is.null(patterns.length))
    patterns.length = max(unlist(lapply(patterns, nrow)), na.rm=TRUE)
  
  res = lapply(patterns, function(q){
    freq = as.numeric(diff(range(index(q))))/(patterns.length-1)
    timeline = seq(min(index(q), na.rm = TRUE), max(index(q), na.rm = TRUE), by=freq)
    res = zoo(data.frame(na.spline(q, xout = timeline)), timeline)
    names(res) = names(q)
    res
  })
  res
}

