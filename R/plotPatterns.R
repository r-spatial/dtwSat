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


#' @title Plotting temporal patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the temporal patterns.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwTimeSeries}}, 
#' \code{\link[zoo]{zoo}}, or list of \code{\link[zoo]{zoo}}.
#' @param labels a vector with labels of the time series. If not declared 
#' the function will plot all time series. 
#' @param year An integer. The base year to shift the dates of the time series to. 
#' If NULL then it does not shif the time series. Default is 2005. 
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting. If not declared the function will plot 
#' all attributes.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}} and 
#' \code{\link[dtwSat]{plotTimeSeries}}
#'  
#' @examples
#' patt = twdtwTimeSeries(patterns.list)
#' plotPatterns(patt)
#' plotPatterns(patt, attr="evi")
#' 
#' @export
plotPatterns = function(x, labels=NULL, attr, year=2005){
  
  if(is(x, "twdtwMatches")) x = x@patterns
  if(is(x, "twdtwTimeSeries")) x = subset(x, labels) 
  x = twdtwTimeSeries(x, labels)
  labels = labels(x)
  
  # Shift dates 
  if(!is.null(year)) x = shiftDates(x, year=year)
    
  # Build data.frame
  if(missing(attr)) attr = names(x[[1]])
  df.p = do.call("rbind", lapply(labels, function(p){
    ts = x[[p]][,attr,drop=FALSE]
    data.frame(Time=index(ts), ts, Pattern=p)
  }))
  df.p = melt(df.p, id.vars=c("Time","Pattern"))
  
  # Plot temporal patterns
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    facet_wrap(~Pattern) + 
    theme(legend.position = "bottom") + 
    scale_x_date(labels = date_format("%b"))
  
  gp
  
}

