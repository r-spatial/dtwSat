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


#' @title Create temporal patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function fits a gam model and retrieves a smoothed 
#' temporal pattern 
#' 
#' @param x A \code{\link[base]{list}} of \code{\link[zoo]{zoo}} such as 
#' retrived by \code{\link[dtwSat]{extractTimeSeries}}. 
#' @param from A character or \code{\link[base]{Dates}} object in the format 
#' "yyyy-mm-dd". If not informed it is equal to the smallest date of the 
#' first element in x 
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} 
#' object in the format "yyyy-mm-dd". If not informed it is equal to the 
#' greatest date of the first element in x 
#' @param attr A vector character or numeric. The attributes in \code{x} to use 
#' @param freq An integer. The frequency of the output patterns 
#' @param formula A formula. Argument to pass to \code{\link[mgcv]{gam}}
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' packege \pkg{mgcv}
#' 
#' 
#' @docType methods
#' 
#' @return A \code{\link[zoo]{zoo}} object 
#' 
#' @examples
#' 
#' ###
#' 
#' 
#' @export
createPattern = function(x, from, to, freq=1, attr, formula, ...){
  
  # Pattern period 
  if( missing(from) | missing(to) ){
    from = as.Date(min(index(x[[1]])))
    to = as.Date(max(index(x[[1]])))
  } 
  from = as.Date(from)
  to = as.Date(to)
  
  # Get formula variables
  if(!is(formula, "formula"))
    stop("missing object formula")
  vars = all.vars(formula)
  
  # Shift dates to match the same period  
  df = do.call("rbind", lapply(x, function(xx){
    res = shiftDates(xx, year=as.numeric(format(to, "%Y")))
    res = window(res, start = from, end = to)
    res = data.frame(time=index(res), res)
  }))
  # names(df)[1] = vars[2]
  
  dates = as.Date(df$time)
  pred_time = seq(from, to, freq)
  
  fun = function(y, ...){
    df = data.frame(y, as.numeric(dates))
    names(df) = vars
    fit = gam(data = df, formula = formula, ...)
    time = data.frame(as.numeric(pred_time))
    names(time) = vars[2]
    predict.gam(fit, newdata = time)
  }
  
  if(missing(attr)) attr = names(df)[-which(names(df) %in% "time")]
  
  res = sapply(as.list(df[attr]), FUN=fun, ...)
  res = zoo(data.frame(res), as.Date(pred_time))
  res
}

