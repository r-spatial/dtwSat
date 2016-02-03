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
#' function to compute the TWDTW local cost matrix [1]
#' 
#' @param alpha numeric. The steepness of logistic weight
#' @param beta numeric. The midpoint of logistic weight 
#' 
#' @docType methods
#' @return An \code{\link[base]{function}} object
#' 
#' 
#' @seealso \code{\link[dtwSat]{twdtw}}
#' 
#' @references 
#' [1] Maus  V,  C\^{a}mara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, X, XX-XX.
#' 
#' @examples
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' log_fun
#' 
#' @export
logisticWeight = function(alpha, beta){
  function(psi) 1 / (1 + exp(1) ^ (alpha * (psi - beta )))
}

