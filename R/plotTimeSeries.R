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


#' @title Plotting time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the temporal patterns.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwTimeSeries}}, 
#' \code{\link[zoo]{zoo}}, or list of \code{\link[zoo]{zoo}}.
#' @param labels a vector with labels of the time series. If missing, all 
#' elements in the list will be plotted (up to a maximum of 6).
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting. If not declared the function will plot 
#' all attributes.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}} and 
#' \code{\link[dtwSat]{plotPatterns}}
#'  
#' @examples
#' ts = twdtwTimeSeries(example_ts.list)
#' plotTimeSeries(ts)
#' plotTimeSeries(ts, attr="evi")
#' 
#' @export
plotTimeSeries = function(x, labels=NULL, attr){
  
  if(is(x, "twdtwMatches")) x = x@timeseries
  if(is(x, "twdtwTimeSeries")) x = subset(x, labels) 
  x = twdtwTimeSeries(x, labels)
  labels = labels(x)
  
  if(length(labels)>6) labels = labels[1:6]
      
  # Build data.frame
  if(missing(attr)) attr = names(x[[1]])
  df.p = do.call("rbind", lapply(labels, function(p){
    ts = x[[p]][,attr,drop=FALSE]
    data.frame(Time=index(ts), ts, Series=p)
  }))
  df.p = melt(df.p, id.vars=c("Time","Series"))
  
  # Plot time series 
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    theme(legend.position = "bottom") + 
    facet_wrap(~Series, scales = "free_x", ncol=1) 
  
  gp
  
}


