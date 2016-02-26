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
#' @param object an twdtwTimeSeries object.
#' @param labels a vector with the patterns labels. If not informed the 
#' function retrieves a pattern for each label in the object twdtwTimeSeries. 
#' 
#' @return an object of class \code{\link[dtwSat]{twdtwTimeSeries}} 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, 
#' \code{\link[dtwSat]{getTimeSeries}}, and 
#' \code{\link[dtwSat]{twdtwApply}}
#' 
#' @export
setGeneric("createPatterns", function(object, ...) standardGeneric("createPatterns"))

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
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' prj_string = scan(system.file("lucc_MT/data/samples_projection", package="dtwSat"), 
#'                   what = "character")
#' 
#' # Extract time series 
#' ts = getTimeSeries(rts, samples = field_samples, proj4string = prj_string)
#' 
#' # Create temporal patterns 
#' createPatterns.twdtwTimeSeries(x=ts, from="2005-09-01", to="2006-09-01", freq=8, formula = y~s(x))
#' 
#' # Plot patterns 
#' autoplot(ts[[1]], facets = NULL) + xlab("Time") + ylab("Value")
#' 
#' 
#' @export
setMethod("createPatterns", 
            function(object, labels, ...) print("hello") )

createPatterns.twdtwTimeSeries = function(x, from=NULL, to=NULL, freq=1, split=TRUE, attr=NULL, formula, ...){

  # Get formula variables
  if(!is(formula, "formula"))
    stop("missing object formula")
  vars = all.vars(formula)
  
  # Split samples according their labels 
  if(split) {
      labels = as.character(labels(x))
      levels = as.character(levels(x))
      names(levels) = levels 
      x = lapply(levels, function(l) x[labels==l] )
  } 
  # Create patterns 
  res = lapply(x, FUN = .createPattern, from=from, to=to, freq=freq, attr=attr, formula=formula, ...)
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
  df = do.call("rbind", lapply(as.list(x), function(x){
    res = shiftDates(x, year=as.numeric(format(to, "%Y")))
    res = res[[1]]
    res = window(res, start = from, end = to)
    res = data.frame(time=index(res), res)
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
  
  if(missing(attr)) attr = names(df)[-which(names(df) %in% "time")]
  
  res = sapply(as.list(df[attr]), FUN=fun, ...)
  labels = unique(labels(x))
  if(length(labels)>1) {
    labels = NULL
    warning("x labels are not unique.")
  }  
  res = twdtwTimeSeries(zoo(data.frame(res), as.Date(pred_time)), labels=labels)
  res
}
