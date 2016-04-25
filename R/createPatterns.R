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
#   R Package dtwSat - 2016-02-18                             #
#                                                             #
###############################################################


#' @title Create patterns 
#' @name createPatterns
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Create temporal patterns from objects of class twdtwTimeSeries.
#' 
#' @param x an object of class \code{\link[dtwSat]{twdtwTimeSeries}}.
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
#' @param split A logical. If TRUE the samples are split by label. If FALSE 
#' all samples are set to the same label. 
#' 
#' @param formula A formula. Argument to pass to \code{\link[mgcv]{gam}}. 
#' 
#' @param ... other arguments to pass to the function \code{\link[mgcv]{gam}} in the 
#' packege \pkg{mgcv}.
#' 
#' @return an object of class \code{\link[dtwSat]{twdtwTimeSeries}} 
#' 
#'
#' @details The hidden assumption is that the temporal pattern is a cycle the repeats itself 
#' within a given time interval. Therefore, all time series samples in \code{x} are aligned 
#' to each other keeping their respective sequence of days of the year. The function fits a 
#' Generalized Additive Model (GAM) to the aligned set of samples.  
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, 
#' \code{\link[dtwSat]{getTimeSeries}}, and 
#' \code{\link[dtwSat]{twdtwApply}}
#' 
#' @export
setGeneric("createPatterns", function(x, ...) standardGeneric("createPatterns"))

#' @rdname createPatterns
#' @aliases createPatterns-twdtwMatches
#' @examples
#' # Creating patterns from objects of class twdtwTimeSeries 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, timeline=timeline)
#' 
#' # Read field samples 
#' \dontrun{
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' prj_string = scan(system.file("lucc_MT/data/samples_projection", package="dtwSat"), 
#'                   what = "character")
#' 
#' # Extract time series 
#' ts = getTimeSeries(rts, y = field_samples, proj4string = prj_string)
#' 
#' # Create temporal patterns 
#' patt = createPatterns(x=ts, from="2005-09-01", to="2006-09-01", freq=8, formula = y~s(x))
#' 
#' # Plot patterns 
#' autoplot(patt[[1]], facets = NULL) + xlab("Time") + ylab("Value")
#' 
#' }
#' @export
setMethod("createPatterns", "twdtwTimeSeries",
          function(x, from=NULL, to=NULL, freq=1, attr=NULL, split=TRUE, formula, ...) {
          
                # Get formula variables
                if(!is(formula, "formula"))
                  stop("missing formula")
                vars = all.vars(formula)
                
                # Split samples according their labels 
                if(split) {
                    levels = as.character(levels(x))
                    labels = as.character(labels(x))
                    names(levels) = levels 
                    x = lapply(levels, function(l) x[labels==l] )
                } else {
                    levels = as.character(levels(x)[1])
                    labels = rep(levels, length(x))
                    names(levels) = levels 
                    x@labels = factor(labels)
                    x = list(x)
                    names(x) = levels
                }
                
                # Create patterns 
                res = lapply(x, FUN = .createPattern, from=from, to=to, freq=freq, attr=attr, formula=formula, ...)
                twdtwTimeSeries(res)
})


.createPattern = function(x, from, to, freq, attr, formula, ...){
  
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
    res = shiftDates(x, year=as.numeric(format(to, "%Y")))
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











