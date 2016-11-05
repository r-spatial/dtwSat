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


#' @title Plotting alignments 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the alignments and TWDTW 
#' dissimilarity measures.
#' 
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwMatches}}.
#' @param timeseries.labels the label or index of the time series.
#' Default is 1. 
#' @param patterns.labels a vector with labels of the patterns. If not 
#' declared the function will plot the alignments for all patterna in \code{x}.
#' @param attr An \link[base]{integer} or \link[base]{character} vector 
#' indicating the attribute for plotting. Default is 1.
#' @param threshold A number. The TWDTW dissimilarity threshold, \emph{i.e.} the 
#' maximum TWDTW cost for consideration. Default is \code{Inf}.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotMatches}}, and
#' \code{\link[dtwSat]{plotClassification}}.
#' 
#' @examples
#' log_fun = logisticWeight(-0.1, 100)
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat1 = twdtwApply(x=ts, y=patt, weight.fun=log_fun)
#' 
#' plotAlignments(mat1)
#' 
#' plotAlignments(mat1, attr=c("evi","ndvi"))
#' 
#' @export
plotAlignments = function(x, timeseries.labels=NULL, patterns.labels=NULL, attr=1, threshold=Inf){
  
  x = subset(x, timeseries.labels[1], patterns.labels) 
    
  ## Get data
  ts = getTimeSeries(x)[[1]]
  alignments = x[[1]]
  
  # Get time series
  df.x = data.frame(ts[,attr,drop=FALSE])
  df.x$Time = as.Date(rownames(df.x))
  df.alignments = melt(df.x, id="Time")
  df.alignments$distance = NA
  df.alignments$variable = as.character(df.alignments$variable)
  df.alignments$Variable = df.alignments$variable
  df.alignments$Pattern = NA
  df.alignments$group = df.alignments$variable 
  df.alignments$facets = 1
  df.alignments$facets = factor(df.alignments$facets, levels = c(1,2), labels = c("Time series","TWDTW dissimilarity measure"))
  
  # Get matching points
  df.matches = list()
  df.matches$Time = c(alignments$from, alignments$to)
  df.matches$variable = rep(alignments$label, 2)
  df.matches$value = rep(alignments$distance, 2)
  df.matches$distance = df.matches$value
  df.matches$Variable = NA
  df.matches$Pattern = df.matches$variable
  df.matches$group = as.character(rep(1:length(alignments$label), 2))
  df.matches$facets = 2
  df.matches = data.frame(df.matches, stringsAsFactors = FALSE)
  df.matches$facets = factor(df.matches$facets, levels = c(1,2), labels = c("Time series","TWDTW dissimilarity measure"))
  
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

