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
#   R Package dtwSat - 2016-01-19                             #
#                                                             #
###############################################################


#' @title Linear weight function 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Builds a linear time weight 
#' function to compute the TWDTW local cost matrix [1].
#' 
#' @param a numeric. The slop of the line.
#' @param b numeric. The intercept of the line. 
#' 
#' @docType methods
#' @return A \code{\link[base]{function}} object.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{\link[dtwSat]{twdtwApply}} through the argument \code{weight.fun}.
#' This will add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful to analyse time series with phenological cycles
#' of vegetation that are usually bound to seasons. In previous studies by 
#' \insertCite{Maus:2016;textual}{dtwSat} the logistic weight had better results than the 
#' linear for land cover classification. See \insertCite{Maus:2016;textual}{dtwSat} and 
#' \insertCite{Maus:2019;textual}{dtwSat}.
#' 
#' @seealso \code{\link[dtwSat]{twdtwApply}}
#' 
#' @references 
#'   \insertAllCited{}
#' 
#' @examples
#' lin_fun = linearWeight(a=0.1)
#' lin_fun
#' 
#' @export
linearWeight = function(a, b=0){
  function(phi, psi) phi + a*psi + b
}

