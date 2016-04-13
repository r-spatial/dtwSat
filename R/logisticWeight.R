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
#' @return An \code{\link[base]{function}} object.
#' 
#' @details The linear \code{linearWeight} and \code{logisticWeight} weight functions 
#' can be passed to \code{\link[dtwSat]{twdtwApply}} through the argument \code{weight.fun}. 
#' This will add a time-weight to the dynamic time warping analysis. The time weight 
#' creates a global constraint useful to analyse time series with phenological cycles
#' of vegetation that are usually bound to seasons. In previous studies by [1] the 
#' logistic weight had better results than the linear for land cover classification. 
#' See [1] for details about the method. 
#' 
#' @seealso \code{\link[dtwSat]{twdtwApply}}
#' 
#' @references 
#' [1] Maus  V,  Camara  G,  Cartaxo  R,  Sanchez  A,  Ramos  FM,  de Queiroz, GR.
#' (2016). A Time-Weighted Dynamic Time Warping method for land use and land cover 
#' mapping. Selected Topics in Applied Earth Observations and Remote Sensing, 
#' IEEE Journal of, vol.PP, no.99, pp.1-11.
#' 
#' @examples
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' log_fun
#' 
#' @export
logisticWeight = function(alpha, beta){
  function(psi) 1 / (1 + exp(alpha * (psi - beta )))
}

