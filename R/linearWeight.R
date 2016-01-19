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
#   R Package dtwSat - 2016-19-01                             #
#                                                             #
###############################################################


#' @title Linear weight function 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Builds a linear time weight 
#' function to compute the TWDTW local cost matrix [1]
#' 
#' @param a numeric. The slop of the line 
#' @param b numeric. The intercept of the line 
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
#' weight.fun = linearWeight(a=0.1)
#' weight.fun
#' 
#' @export
linearWeight = function(a, b=0){
  function(psi) a*psi + b
}

