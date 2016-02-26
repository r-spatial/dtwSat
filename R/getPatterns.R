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
#   R Package dtwSat - 2016-02-18                             #
#                                                             #
###############################################################

#' @export
setGeneric("getPatterns", function(object, labels=NULL) standardGeneric("getPatterns"))

#' @inheritParams twdtwMatches-class
#' @describeIn twdtwMatches Get patterns from objects of class twdtwMatches.
#'
#' @examples 
#' # Getting patterns from objects of class twdtwMatches
#' patt = twdtwTimeSeries(patterns.list)
#' ts = twdtwTimeSeries(example_ts.list)
#' mat = twdtwApply(x=ts, y=patt, weight.fun=logisticWeight(-0.1,50))
#' getPatterns(mat)
#' getTimeSeries(mat)
#' 
#' @export
setMethod("getPatterns", c("twdtwMatches","ANY"),
          function(object, labels) getPatterns.twdtwMatches(object, labels) )

getPatterns.twdtwMatches = function(object, labels) {
  getTimeSeries(object@patterns, labels)
}

