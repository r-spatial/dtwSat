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


#' @title Logistic weight function 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Builds a logistic time weight 
#' function to compute the TWDTW local cost matrix
#' 
#' @param alpha numeric. The steepness of logistic weight
#' @param beta numeric. The midpoint of logistic weight
#' @param theta numeric between 0 and 1. The weight of the time 
#' for the TWDTW computation. Use \code{theta=0} to cancel the time-weight, 
#' \emph{i.e.} to run the original DTW algorithm. Default is 0.5. 
#' For details see [1]
#' 
#' @docType methods
#' @return An \code{\link[base]{function}} object
#' 
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}
#' 
#' @references 
#' [1] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2015). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, X, XX-XX.
#' 
#' @examples
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' weight.fun
#' 
#' @export
logisticWeight = function(alpha, beta, theta=0.5){
  function(phi, psi) (1-theta)*phi + theta*(1 / (1 + exp(1) ^ (alpha * (psi - beta ))))
}

