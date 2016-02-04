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

#' @title Get pattern names from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves pattern names in the 
#' \link[dtwSat]{twdtw-class} object.
#' 
#' @param object A \link[dtwSat]{twdtw-class} object.
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves the names for all patterns.
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{character}}
#' or \code{\link[base]{numeric}} vector.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}.
#' 
#' @examples
#' 
#' matches = twdtw(x=example_ts, patterns=patterns.list)
#' 
#' getPatternNames(matches)
#' 
#' getPatternNames(matches, p.names=c(1,3))
#' 
#' getPatternNames(matches, p.names="Maize")
#' 
#' @export
setGeneric("getPatternNames", 
           function(object, p.names){
             p.names = .getPatternNames(object, p.names)
             if(any(is.na(p.names)))
               warning("the patterns identification is invalid", call. = FALSE)
             p.names
           }
)

.getPatternNames = function(object, p.names){
  if(missing(p.names)) p.names = seq_along(object@alignments)
  all_names = names(object@alignments)
  names(all_names) = all_names
  if(is.null(all_names)) all_names = seq_along(object@alignments)
  all_names[p.names]
}

