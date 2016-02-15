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


#' @title Plotting subintervals classification  
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the classification of each 
#' subinterval of the time series based on TWDTW analysis. 
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object.
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting. If not declared the function will plot 
#' all attributes.
#' @param ... additional arguments passed to \code{\link[dtwSat]{classifyIntervals}}.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotMatches}}, and 
#' \code{\link[dtwSat]{plotPatterns}}.
#' 
#' @examples
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun, 
#'         normalize.patterns=TRUE, patterns.length=23)
#' 
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' by = "6 month"
#' 
#' # All classes
#' gp = plotClassification(x=matches, from=from, to=to, by=by, overlap=.4)
#' gp
#' 
#' 
#' @export
plotClassification = function(x, attr, ...){
  
  ## Get data
  ts = getTimeSeries(x)
  best_class = classifyIntervals(x, simplify = FALSE, ...)
  tx = index(ts)
  
  I = min(best_class$from, na.rm = TRUE)-30 <= index(ts) & 
    index(ts) <= max(best_class$to, na.rm = TRUE)+30
  
  if(missing(attr)) attr=names(ts)
  xx = ts[I,attr,drop=FALSE]
  tx = index(xx)
  
  df.x = melt(data.frame(Time=tx, xx), id="Time")
  
  y.labels = pretty_breaks()(range(df.x$value, na.rm = TRUE))
  y.breaks = y.labels
  df.pol = do.call("rbind", lapply(1:nrow(best_class), function(i){
    data.frame(
      Time = c(best_class$from[i], best_class$to[i], best_class$to[i], best_class$from[i]),
      Group = rep(i, 4),
      Class = rep(best_class$pattern[i], 4),
      value = rep(range(y.breaks, na.rm = TRUE), each=2))
  }))
  df.pol$Group = factor(df.pol$Group)
  df.pol$Class = factor(df.pol$Class)
  
  gp = ggplot() +
    geom_polygon(data=df.pol, aes_string(x='Time', y='value', 
                                         group='Group', fill='Class'), alpha=.7) +
    scale_fill_brewer(palette="Set3") + 
    geom_line(data=df.x, aes_string(x='Time', y='value', colour='variable')) +
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab("Value") + 
    xlab("Time")
  gp
}

