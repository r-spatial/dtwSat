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
#   R Package dtwSat - 2016-08-12                             #
#                                                             #
###############################################################


#' @title class "twdtwAccuracy" 
#' @name twdtwAccuracy-class
#' @aliases twdtwAccuracy
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This class stores the cross-validation. 
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwRaster}}.
#' 
#' @param x an object of class \code{\link[base]{data.frame}}.
#' 
#' @seealso 
#' \code{\link[dtwSat]{createPatterns}}, and 
#' \code{\link[dtwSat]{twdtwApply}}.
#'
#' @section Slots :
#' \describe{
#'  \item{\code{samples}:}{A data frame with mapped and reference classes.}
#'  \item{\code{accuracy}:}{A list with the accuracy metrics.}
#' }
#'
#' @examples
#' \dontrun{
#' 
#' }
NULL

df = do.call("rbind", lapply(twdtw_res[], function(xx) xx[which.min(xx$distance),]) )

data = data.frame(.adjustFactores(ref, pred), df[,!names(df)%in%"labels"])
# error.matrix = table(Predicted=data$Predicted, Reference=data$Reference)



