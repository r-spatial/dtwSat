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


#' @title Plotting matching points 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the matching points from   
#' TWDTW analysis.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwMatches}}.
#' @param timeseries.labels the label or index of the time series.
#' Default is 1. 
#' @param patterns.labels a vector with labels of the patterns. If not 
#' declared the function will plot one alignment for each pattern.
#' @param k A positive integer. The index of the last alignment to include in 
#' the plot. If not declared the function will plot the best match for 
#' each pattern. 
#' @param attr An \link[base]{integer} or \link[base]{character} vector 
#' indicating the attribute for plotting. Default is 1.
#' @param shift A number, it shifts the pattern position in the \code{x}
#' direction. Default is 0.5.
#' @param show.dist show the distance for each alignment. Default is FALSE.
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, and
#' \code{\link[dtwSat]{plotClassification}}.
#'   
#' @examples
#' log_fun = logisticWeight(-0.1, 100)
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat1 = twdtwApply(x=ts, y=patt, weight.fun=log_fun, keep=TRUE)
#' 
#' plotMatches(mat1)
#' 
#' plotMatches(mat1, patterns.labels="Soybean", k=4)
#' 
#' plotMatches(mat1, patterns.labels=c("Soybean","Maize"), k=4)
#' 
#' plotMatches(mat1, patterns.labels=c("Soybean","Cotton"), k=c(3,1))
#' 
#' @export
plotMatches = function(x, timeseries.labels=1, patterns.labels=NULL, k=1, attr=1, shift=0.5, show.dist=FALSE){
 
  x = subset(x, timeseries.labels[1], patterns.labels) 
  ## Get data
  internals = getInternals(x)[[1]]
  if(any(sapply(internals, function(x) length(x$internals))<1))
    stop("plot methods requires internals, set keep=TRUE on twdtwApply() call")
  matching = getMatches(x)[[1]]
  alignments = getAlignments(x)[[1]]
  ts = getTimeSeries(x)[[1]]
  patterns = getPatterns(x)
  
  y = as.character(labels(x)$patterns)
  if(length(k)==1){
    y = rep(y, each=k)
    k = unlist(lapply(table(y), function(i) seq(from=1, to=i) ))
  }
  if(length(y)!=length(k))
    stop("if length of k greater than 1, then patterns.labels must have the same length as k.")
  
  xx = ts[,attr,drop=FALSE]
  tx = index(xx)
  
  y.labels = pretty_breaks()(range(xx, na.rm = TRUE))
  y.breaks = y.labels
  
  # Get time series 
  df.x = data.frame(Time=tx, xx)
  
  # Build matching points data.frame
  df.list = lapply(seq_along(y), function(i){
    p = y[i]
    yy = patterns[[p]][,attr,drop=FALSE]
    ty = index(yy)
    
    if(k[i]>alignments[[p]]$K){
      warning("alignment index out of bounds", call. = TRUE)
      return(NULL)
    }
      
    map = data.frame(matching[[p]]$matching[[k[i]]])
    delay = tx[map$index2[1]]-ty[1]
    if(delay>0)
      delay = delay + diff(range(ty, na.rm = TRUE))*shift
    if(delay<0)
      delay = delay - diff(range(ty, na.rm = TRUE))*shift
    
    df.pt = data.frame(Time=ty[map$index1]+delay, yy[map$index1,,drop=FALSE]+max(xx, na.rm = TRUE))
    df.match.pt = df.pt
    df.match.pt$alig = paste(1:nrow(map),p,k[i],sep="_")
    df.match.x = df.x[map$index2,]
    df.match.x$alig = paste(1:nrow(map),p,k[i],sep="_")
    df.match = rbind(df.match.pt, df.match.x)
    df.pt$Matches = paste(p,k[i])
    df.dist = data.frame(Time=max(ty[map$index1]+delay)+diff(range(df.pt$Time))/3,
                         max(df.pt[,names(yy)]),Dist=alignments[[p]]$distance[k[i]])
    names(df.dist) = c("Time", names(yy), "Dist")
    list(match=df.match, pt=df.pt, dist=df.dist)
  })
  
  df.pt = do.call("rbind", lapply(df.list, function(df) df$pt))
  df.match = do.call("rbind", lapply(df.list, function(df) df$match))
  
  attr_names = names(df.x)[2]
  gp = ggplot(data=df.x, aes_string(x='Time', y=eval(attr_names))) +
    geom_line() +
    geom_line(data=df.pt, aes_string(x='Time', y=eval(attr_names), 
                                     group='Matches', colour='Matches')) + 
    geom_line(data=df.match, linetype = 2, colour = "grey", 
              aes_string(x='Time', y=eval(attr_names), group='alig')) + 
    scale_y_continuous(breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(attr_names) 
  
  if(show.dist){
    df.dist = do.call("rbind", lapply(df.list, function(df) df$dist))
    df.dist$Dist = paste("Distance:",round(df.dist$Dist,2))
    gp = gp + geom_text(data=df.dist, mapping = aes_string(x='Time', y=eval(attr_names), label='Dist')) 
  }
  
  gp
}


