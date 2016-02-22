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
#   R Package dtwSat - 2016-01-19                             #
#                                                             #
###############################################################


#' @title Joint TWDTW alignments 
#' @name joinAlignments
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Join TWDTW alignments from two or more objects of class twdtwMatches. 
#' 
#' @param ... objects of class twdtwMatches 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @examples 
#' ## 
#' 
#' @export joinAlignments 
setGeneric(name = "joinAlignments",  
          def = function(...) standardGeneric("joinAlignments")
)

#' @inheritParams joinAlignments
#' @describeIn twdtwMatches Join TWDTW alignments from two or more objects of class twdtwMatches 
setMethod("joinAlignments", "twdtwMatches",
          definition = function(...) joinAlignments.twdtwMatches(list(...)))
         
joinAlignments.twdtwMatches = function(x){
  n = length(x)
  if(n==1) return(x[[1]])
  align = x[[n]]@alignments
  if(length(align)<1) return(joinAlignments.twdtwMatches(x[-n]))
  x[[1]]@alignments = c(x[[1]]@alignments, align)
  x[[1]]@timeseries = joinTimeSeries(x[[1]]@timeseries, x[[n]]@timeseries)
  x[[1]]@patterns = joinPatterns(x[[1]]@patterns, x[[n]]@patterns)
  joinAlignments.twdtwMatches(x[-n])
}



          