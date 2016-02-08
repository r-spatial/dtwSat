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

#' @title Show method for twdtw-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Show the head of the matches ordered by 
#' TWDTW similarity.
#' 
#' @param object An \code{\link[dtwSat]{twdtw-class}} object.
#' 
#' @docType methods
#' 
#' @return Print the twdtw-class object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and 
#' \code{\link[dtwSat]{twdtw}}.
#'  
#' @examples
#' 
#' matches = twdtw(x=example_ts, patterns=patterns.list)
#'        
#' show(matches)
#' 
#' @rdname show-method
#' 
#' @export
setMethod("show", 
          signature = signature(object="twdtw"),
          definition = function(object){
            cat("Time-Weighted DTW alignment object\n")
            cat("Number of alignments:",nrow(getAlignments(object)),"\n")
            print(head(getAlignments(object)))
            invisible(NULL)
          }
)

