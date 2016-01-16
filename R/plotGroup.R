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


#' @title Plotting time interval groups 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the group of each 
#' time intervals based on \code{\link[dtwSat]{twdtw}} analysis 
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting, \emph{i.e.} a column of the \code{x}. 
#' If not declared the function will plot all attributes
#' @param ... additional arguments passed to \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @return A \link[ggplot2]{ggplot} object
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @examples
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' # Classify interval
#' from = as.Date("2009-09-01")
#' to = as.Date("2013-09-01")
#' by = "6 month"
#' 
#' # All classes
#' gp = plotGroup(x=alig, from=from, to=to, by=by, overlap=.3)
#' gp
#' 
#' # Cotton and Maize 
#' gp = plotGroup(x=alig, attr=c("ndvi","evi"), from=from, to=to, by=by, 
#'              overlap=.3, p.names=c("Cotton","Maize"))
#' gp
#' 
#' 
#' @export
plotGroup = function(x, attr, ...){
  
  ## Get data
  internals = getInternals(x, 1)[[1]]
  if(is.null(internals))
    stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
  
  best_class = classifyIntervals(x, ...)
  # best_class = classifyIntervals(x, from=from, to=to, by=by, overlap=.3, threshold=Inf)
  
  tx = index(internals$x)
  
  I = min(best_class$from, na.rm = TRUE)-30 <= index(internals$x) & 
    index(internals$x) <= max(best_class$to, na.rm = TRUE)+30
  
  if(missing(attr)) attr=names(internals$x)
  xx = internals$x[I,attr,drop=FALSE]
  tx = index(xx)
  
  df.x = melt(data.frame(Time=tx, xx), id="Time")
  
  y.labels = pretty_breaks()(range(df.x$value, na.rm = TRUE))
  y.breaks = y.labels
  df.pol = do.call("rbind", lapply(1:nrow(best_class), function(i){
    data.frame(
      Time = c(best_class$from[i], best_class$to[i], best_class$to[i], best_class$from[i]),
      Group = rep(i, 4),
      Pattern = rep(best_class$pattern[i], 4),
      value = rep(range(y.breaks, na.rm = TRUE), each=2))
  }))
  df.pol$Group = factor(df.pol$Group)
  df.pol$Pattern = factor(df.pol$Pattern)
  
  gp = ggplot() +
    geom_polygon(data=df.pol, aes_string(x='Time', y='value', 
                                         group='Group', fill='Pattern'), alpha=.7) +
    scale_fill_brewer(palette="Set3") + 
    geom_line(data=df.x, aes_string(x='Time', y='value', colour='variable')) +
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab("Value") + 
    xlab("Time")
  gp
}

