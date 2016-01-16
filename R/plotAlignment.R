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


#' @title Plotting twdtw alignments
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the alignments and respective TWDTW 
#' distances 
#' 
#' 
#' @param x An \code{\link[dtwSat]{twdtw-class}} object
#' @param p.names A \link[base]{character} or \link[base]{integer}
#' vector with the patterns identification. If not declared the function 
#' will plot alignments for all patterna in the 
#' \code{\link[dtwSat]{twdtw-class}} object
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting, \emph{i.e.} a column of the \code{pattern}. 
#' Default is 1
#' @param threshold A number. The TWDTW threshold, i.e. the maximum TWDTW 
#' cost for consideration. Default is \code{Inf}
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' gp = plotAlignment(alig, attr=c("ndvi","evi"))
#' gp
#' 
#' gp = plotAlignment(alig, attr=c("ndvi","evi"), threshold=4)
#' gp
#' 
#' gp = plotAlignment(x=alig, p.names="Cotton", 
#'                    attr=c("ndvi","evi"), threshold=6)
#' gp
#' 
#' @export
plotAlignment = function(x, p.names, attr=1, threshold=Inf){
  
  if(missing(p.names)) {
    p.names = getPatternNames(x)
  } else {
    p.names = getPatternNames(x, p.names)
  }
  
  ## Get data
  internals = getInternals(x, 1)[[1]]
  if(is.null(internals))
    stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
  
  alignments = getAlignments(x, p.names)
  
  # Get time series
  df.x = data.frame(internals$x[,attr,drop=FALSE])
  df.x$Time = as.Date(rownames(df.x))
  df.alignments = melt(df.x, id="Time")
  df.alignments$distance = NA
  df.alignments$variable = as.character(df.alignments$variable)
  df.alignments$Variable = df.alignments$variable
  df.alignments$Pattern = NA
  df.alignments$group = df.alignments$variable 
  df.alignments$facets = 1
  df.alignments$facets = factor(df.alignments$facets, levels = c(1,2), labels = c("Time series","TWDTW alignments distance"))
  
  # Get matching points
  df.matches = list()
  df.matches$Time = c(alignments$from, alignments$to)
  df.matches$variable = rep(alignments$pattern, 2)
  df.matches$value = rep(alignments$distance, 2)
  df.matches$distance = df.matches$value
  df.matches$Variable = NA
  df.matches$Pattern = df.matches$variable
  df.matches$group = as.character(rep(1:length(alignments$pattern), 2))
  df.matches$facets = 2
  df.matches = data.frame(df.matches, stringsAsFactors = FALSE)
  df.matches$facets = factor(df.matches$facets, levels = c(1,2), labels = c("Time series","TWDTW alignments distance"))
  
  I = which(df.matches$value>threshold)
  if(length(I)>0)
    df.matches = df.matches[-I,]
  
  df.all = rbind(df.alignments, df.matches)  
  
  gp = ggplot(data=df.all) +  
    geom_line(data=df.alignments, aes_string(x='Time', y='value', group='group', linetype='Variable')) + 
    facet_wrap(~facets, ncol = 1, scales = "free_y") + 
    geom_path(data=df.matches, aes_string(x='Time', y='distance', group='group', colour='Pattern')) +
    geom_point(data=df.matches, aes_string(x='Time', y='distance', group='group', colour='Pattern')) 
  gp
}

