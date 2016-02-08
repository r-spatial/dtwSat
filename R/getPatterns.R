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
#   R Package dtwSat - 2016-02-08                             #
#                                                             #
###############################################################

#' @title Get patterns from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves temporal patterns used in the
#' TWDTW analysis.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves all patterns. 
#' 
#' @docType methods
#' @return A a list of \link[zoo]{zoo} objects
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}.
#' 
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, 
#'              weight.fun = log_fun)
#' 
#' a = getPatterns(matches)
#' names(a) 
#' 
#' a = getPatterns(matches, p.names="Maize")
#' names(a) 
#' 
#' a = getPatterns(matches, p.names=c(1,2))
#' names(a) 
#' 
#' 
#' @export
setGeneric("getPatterns", 
           function(object, p.names){
             p.names = getPatternNames(object, p.names)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             object@patterns[p.names]
           }
)

