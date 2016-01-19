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

#' @title Get alignments from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the alignments 
#' from the object \link[dtwSat]{twdtw-class}
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \link[base]{data.frame} with the following columns:  
#'       \cr\code{pattern}: the pattern identification,
#'       \cr\code{from}: starting date,
#'       \cr\code{to}: ending date, and
#' 	     \cr\code{distance}: TWDTW distances.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun)
#' 
#' getAlignments(alig)
#' 
#' getAlignments(alig, p.names="Soybean")
#' 
#' getAlignments(alig, p.names=c(2,3))
#' 
#' @export
setGeneric("getAlignments", 
           function(object, ...) {
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getAlignments(object, p.names)
           }
)

.getAlignments = function(object, p.names){
  k = 1
  names(p.names) = NULL
  do.call("rbind", lapply(p.names, function(p){
    x = object@alignments[[p]]
    name = numeric(0)
    r.names = NULL
    if(length(x$distance)>0){
      name = p
      r.names = paste0(k:(k+x$K-1))
      k <<- k + x$K
    }
    data.frame(pattern  = name,
               from     = x$from,
               to       = x$to,
               distance = x$distance,
               stringsAsFactors = FALSE,
               row.names = r.names
    )
  }))
}

