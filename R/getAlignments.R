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

#' @title Get alignments from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the TWDTW alignments.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' @param y A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves the alignments for all patterns. 
#' 
#' @docType methods
#' @return A \link[base]{data.frame} with the following attributes:  
#'       \cr\code{pattern}: the pattern identification,
#'       \cr\code{from}: starting date,
#'       \cr\code{to}: ending date, and
#' 	     \cr\code{distance}: TWDTW dissimilarity.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}.
#' 
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun)
#' 
#' getAlignments(matches)
#' 
#' getAlignments(matches, y="Soybean")
#' 
#' getAlignments(matches, y=c(2,3))
#' 
#' @export
setGeneric("getAlignments", 
           function(object, y) {
             y = getPatternNames(object, y)
             if(any(is.na(y)))
               stop("the patterns identification is invalid")
             .getAlignments(object, y)
           }
)

.getAlignments = function(object, y){
  k = 1
  names(y) = NULL
  do.call("rbind", lapply(y, function(p){
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

