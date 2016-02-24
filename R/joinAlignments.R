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
#' @return an object of class \code{\link[dtwSat]{twdtwMatches}}
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, and 
#' \code{\link[dtwSat]{joinTimeSeries}}
#' 
#' @export  
setGeneric(name = "joinAlignments",  
          def = function(...) standardGeneric("joinAlignments")
)


#' @rdname joinAlignments
#' @aliases joinAlignments-twdtwMatches
#' @examples
#' # Join TWDTW align from objects of class twdtwTimeSeries
#' patterns = twdtwTimeSeries(timeseries=patterns.list, labels=names(patterns.list))
#' ts1 = twdtwTimeSeries(timeseries=example_ts.list[1], labels="A")
#' ts2 = twdtwTimeSeries(timeseries=example_ts.list[2], labels="B")
#' mat1 = twdtwApply(x=ts1, y=patterns)
#' mat2 = twdtwApply(x=ts2, y=patterns)
#' matches = joinAlignments(mat1, mat2)
#' matches
#' 
#' @export
setMethod("joinAlignments", "twdtwMatches",
          definition = function(...) joinAlignments.twdtwMatches(list(...)))

joinAlignments.twdtwMatches = function(x){
  n = length(x)
  if(n==1) return(x[[1]])
  align = x[[n]]@alignments
  if(length(align)<1) return(joinAlignments.twdtwMatches(x[-n]))
  x[[1]]@alignments = c(x[[1]]@alignments, align)
  x[[1]]@timeseries = joinTimeSeries(x[[1]]@timeseries, x[[n]]@timeseries, join.labels=FALSE)
  x[[1]]@patterns = joinTimeSeries(x[[1]]@patterns, x[[n]]@patterns, join.labels=TRUE)
  joinAlignments.twdtwMatches(x[-n])
}



          