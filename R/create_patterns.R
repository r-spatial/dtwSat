#' @title Create patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Create temporal patterns from objects of class twdtwTimeSeries.
#' 
#' @param x an object of class \code{\link[base]{data.frame}}.
#' 
#' @param from A character or \code{\link[base]{Dates}} object in the format 
#' "yyyy-mm-dd". If not provided it is equal to the smallest date of the 
#' first element in x. See details. 
#'  
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} 
#' object in the format "yyyy-mm-dd". If not provided it is equal to the 
#' greatest date of the first element in x. See details. 
#' 
#' @param attr A vector character or numeric. The attributes in \code{x} to be used. 
#' If not declared the function uses all attributes.  
#' 
#' @param freq An integer. The sampling frequency of the output patterns.
#'
#' @param split A logical. If TRUE the samples are split by label. If FALSE 
#' all samples are set to the same label. 
#' 
#' @param formula A formula. Argument to pass to \code{\link[mgcv]{gam}}. 
#' 
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' package \pkg{mgcv}.
#' 
#' @return an object of class \code{\link[base]{data.frame}} 
#' 
#'
#' @details The hidden assumption is that the temporal pattern is a cycle the repeats itself 
#' within a given time interval. Therefore, all time series samples in \code{x} are aligned 
#' with each other, keeping their respective sequence of days of the year. The function fits a 
#' Generalized Additive Model (GAM) to the aligned set of samples.  
#' 
#' @export
create_pattern = function(x, from, to, freq, attr, formula, ...){
  
  # Pattern period 
  if( is.null(from) | is.null(to) ){
    from = as.Date(min(index(x[[1]])))
    to = as.Date(max(index(x[[1]])))
  }

  from = as.Date(from)
  to = as.Date(to)
  
  # Get formula variables
  vars = all.vars(formula)
  
  # Shift dates to match the same period  
  df = do.call("rbind", lapply(as.list(x), function(x){
    res = shift_dates(x, year=as.numeric(format(to, "%Y")))
    res = window(res, start = from, end = to)
    res = data.frame(time=index(res), res)
    names(res) = c("time", names(x))
    res
  }))
  names(df)[1] = vars[2]
  
  dates = as.Date(df[[vars[2]]])
  pred_time = seq(from, to, freq)
  
  fun = function(y, ...){
    df = data.frame(y, as.numeric(dates))
    names(df) = vars
    fit = gam(data = df, formula = formula, ...)
    time = data.frame(as.numeric(pred_time))
    names(time) = vars[2]
    predict.gam(fit, newdata = time)
  }
  
  if(is.null(attr)) attr = names(df)[-which(names(df) %in% vars[2])]

  res = sapply(as.list(df[attr]), FUN=fun, ...)
  zoo(data.frame(res), as.Date(pred_time))
}











