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
#   R Package dtwSat - 2016-30-01                             #
#                                                             #
###############################################################

#' @title Get the number of alignments from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the number of alignments 
#' from the object \link[dtwSat]{twdtw-class}
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \link[base]{numeric} 
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
#' nalignments(alig)
#' 
#' nalignments(alig, p.names="Soybean")
#' 
#' nalignments(alig, p.names=c(2,3))
#' 
#' @export
setGeneric("nalignments", 
           function(object, ...) {
               x = getAlignments(object, ...)
               nrow(x)
           }
)


