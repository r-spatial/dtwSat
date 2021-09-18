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


#' @title Logistic weight function 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Builds a logistic time weight 
#' function to compute the TWDTW local cost matrix [1].
#' 
#' @param alpha numeric. The steepness of logistic weight.
#' @param beta numeric. The midpoint of logistic weight.
#' 
#' @docType methods
#' @return A \code{\link[base]{function}} object.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{\link[dtwSat]{twdtwApply}} through the argument \code{weight.fun}. 
#' This will add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful to analyze time series with phenological cycles
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
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' log_fun
#' 
#' @export
logisticWeight = function(alpha, beta){
  function(phi, psi) phi + 1 / (1 + exp(alpha * (psi - beta )))
}

