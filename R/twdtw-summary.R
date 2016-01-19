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

#' @title Summary method for twdtw-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Summary the alignments for each pattern in the 
#' twdtw-class object 
#' 
#' @param object An \code{\link[dtwSat]{twdtw-class}} object
#' @param ... additional arguments passed to summary
#' 
#' @docType methods
#' 
#' @return A \link[base]{data.frame} object with the the summary 
#' for each pattern 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and 
#' \code{\link[dtwSat]{twdtw}}
#'  
#' @examples
#' 
#' alig = twdtw(x=template, patterns=patterns.list)
#'        
#' show(alig)
#' summary(alig)
#' 
#' @export
setMethod("summary", 
          signature(object = "twdtw"),
          function(object, ...){
            summary.twdtw(object, ...)
          }
)

summary.twdtw = function(object, ...){
  res1 = do.call("rbind", lapply(object@alignments, function(pattern){
    c(N.Alig=length(pattern$distance), summary(pattern$distance))
  }))
  data.frame(res1)
}

