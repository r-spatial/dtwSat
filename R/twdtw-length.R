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

#' @title Length method for twdtw-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Retrieves the number of matches in the 
#' twdtw-class object.
#' 
#' @param x An \code{\link[dtwSat]{twdtw-class}} object.
#' 
#' @docType methods
#' 
#' @return An \link[base]{integer}.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and 
#' \code{\link[dtwSat]{twdtw}}
#'  
#' @examples
#' 
#' matches = twdtw(x=example_ts, patterns=patterns.list)
#'        
#' length(matches)
#' 
#' @rdname length-method
#' @name length
#' 
#' @export
length.twdtw = function(x){
  nrow(getAlignments(x))
}

