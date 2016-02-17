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


#' @title Create temporal patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function fits a Generalized Additive Model model 
#' to the samples and retrieves smoothed temporal patterns. 
#' 
#' @param x A \code{\link[base]{list}} of \code{\link[zoo]{zoo}} such as 
#' retrived by \code{\link[dtwSat]{extractTimeSeries}}. 
#' 
#' @param ref A vector identifying each time series. It must have the same length 
#' as \code{x}. If \code{ref} is not declared then all samples are taken as belonging 
#' to the same class. 
#' 
#' @param from A character or \code{\link[base]{Dates}} object in the format 
#' "yyyy-mm-dd". If not informed it is equal to the smallest date of the 
#' first element in x. See details. 
#'  
#' @param to A \code{\link[base]{character}} or \code{\link[base]{Dates}} 
#' object in the format "yyyy-mm-dd". If not informed it is equal to the 
#' greatest date of the first element in x. See details. 
#' 
#' @param attr A vector character or numeric. The attributes in \code{x} to be used. 
#' If not declared the function uses all attributes.  
#' 
#' @param freq An integer. The sampling frequency of the output patterns.
#'  
#' @param formula A formula. Argument to pass to \code{\link[mgcv]{gam}}. 
#' 
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' packege \pkg{mgcv}.
#' 
#' @details The hidden assumption is that the temporal pattern is a cycle the repeats itself 
#' within a given time interval. Therefore, all time series samples in \code{x} are aligned 
#' to each other keeping the time sequence of each sample and its respective day of the year. 
#' After that the function fits a Generalized Additive Model to the samples.  
#' 
#' @docType methods
#' 
#' @return A \code{\link[zoo]{zoo}} object 
#' 
#' @examples
#' 
#' \dontrun{
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package = "dtwSat"))
#' load(system.file("lucc_MT/field_samples_ts.RData", package = "dtwSat"))
#' 
#' # Sellect Soybean-cotton samples  
#' I = field_samples$class=="Soybean-cotton"
#' soybean_cotton_samples = field_samples_ts[I]
#' 
#' p_soybean_cotton = createPatterns(x = soybean_cotton_samples, 
#'                    formula = y ~ s(x))
#' 
#' plotPatterns(p_soybean_cotton)
#' 
#' p = createPatterns(x = field_samples_ts, ref = field_samples$class,
#'                    formula = y ~ s(x))
#' plotPatterns(p)
#' 
#' }
#' @export
createPatterns = function(x, ref, from, to, freq=1, attr, formula, ...){

  # Get formula variables
  if(!is(formula, "formula"))
    stop("missing object formula")
  vars = all.vars(formula)
  
  if(missing(ref)) ref = rep("class", length(x))
  if(length(ref)!=length(x)) 
    stop("ref and x must have the same length")
    
  unique_ref = as.character(unique(ref))
  names(unique_ref) = unique_ref  
  
  
  # Split samples time series according to the class 
  xx = lapply(unique_ref, function(r) x[ref==r] )

  # Apply createPattern for each unique reference  
  res = lapply(xx, FUN = .createPattern, from=from, to=to, freq=freq, attr=attr, formula=formula, ...)
  res
}


.createPattern = function(x, from, to, freq, attr, formula, ...){
  
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
    res
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

