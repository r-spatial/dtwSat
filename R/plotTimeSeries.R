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
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}} and 
#' \code{\link[dtwSat]{plotPatterns}}
#'  
#' @examples
#' ts = twdtwTimeSeries(example_ts)
#' plotTimeSeries(ts)
#' 
#' @export
plotTimeSeries = function(x, labels=NULL){
  
  if(is(x, "twdtwTimeSeries")) x = subset(x, labels) 
  x = twdtwTimeSeries(x, labels)
  labels = labels(x)
  
  if(length(labels)>6) labels = labels[1:6]
      
  # Build data.frame
  df.p = do.call("rbind", lapply(labels, function(p)
    data.frame(Time=index(x[[p]]), x[[p]], Series=p)
  ))
  df.p = melt(df.p, id.vars=c("Time","Series"))
  
  # Plot time series 
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    facet_wrap(~Series, scales = "free_x") 
  
  gp
  
}


