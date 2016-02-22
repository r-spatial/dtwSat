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
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @include getTimeSeries.R
#' @title Get TWDTW alignments 
#' @name getAlignments
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Generic method to get the TWDTW alignments from objects of class twdtwMatches 
#' 
#' @inheritParams twdtwTimeSeries-class
#' @param ... other arguments to be passed to methods 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and
#' \code{\link[dtwSat]{twdtwRaster-class}}
#' 
#' @examples 
#' matches = new("twdtwMatches")
#' getAlignments(matches)
#' 
#' @export
setGeneric("getAlignments", function(object) standardGeneric("getAlignments"))


#' @inheritParams getAlignments
#' @describeIn twdtwMatches Get the TWDTW alignments from objects of class twdtwMatches 
setMethod(f = "getAlignments", signature = "twdtwMatches", 
          definition = function(object) getAlignments.twdtwMatches(object) )


getAlignments.twdtwMatches = function(object){
  if(length(object@alignments)<1) return(list())
  res = lapply(object@alignments, function(xx){
    if(length(xx)<1) return(NULL)
    res = do.call("rbind", lapply(seq_along(xx), function(j){
      x = xx[j]
      if(x[[j]]$K<1) return(NULL)
      data.frame(pattern  = names(x[j]),
                 from     = x[[j]]$from,
                 to       = x[[j]]$to,
                 distance = x[[j]]$distance,
                 stringsAsFactors = FALSE
      )
    }))
    row.names(res) = as.character(1:nrow(res))
  })
  res
}

