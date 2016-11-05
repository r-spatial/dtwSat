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


#' @title Plotting paths 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting low cost paths in the TWDTW 
#' cost matrix.
#' 
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwMatches}}.
#' @param timeseries.labels the label or index of the time series.
#' Default is 1. 
#' @param patterns.labels a vector with labels of the patterns. If not 
#' declared the function will plot one alignment for each pattern.
#' @param matrix.name A character. The name of the matrix to plot,
#' "costMatrix" for accumulated cost, "localMatrix" for local cost, 
#' or "timeWeight" for time-weight. Default is "costMatrix".
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}}, 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotPaths}},
#' \code{\link[dtwSat]{plotMatches}}, and
#' \code{\link[dtwSat]{plotClassification}}.
#' 
#' @examples
#' log_fun = logisticWeight(-0.1, 100)
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat1 = twdtwApply(x=ts, y=patt, weight.fun=log_fun, keep=TRUE)
#' 
#' plotCostMatrix(mat1, matrix.name="costMatrix")
#' 
#' plotCostMatrix(mat1, matrix.name="localMatrix")
#' 
#' plotCostMatrix(mat1, matrix.name="timeWeight")
#' 
#' @export
plotCostMatrix = function(x, timeseries.labels=NULL, patterns.labels=NULL, matrix.name="costMatrix"){
  
  pt = pmatch(matrix.name,c("costMatrix", "localMatrix", "timeWeight"))
  if(is.na(pt))
    stop("matrix.name is not costMatrix, localMatrix, or timeWeight")
  
  legend_name = c("Accumulated cost", "Local cost", "Time weight")[pt]
  
  x = subset(x, timeseries.labels[1], patterns.labels) 
  y = as.character(labels(x)$patterns)
  ## Get data
  internals = getInternals(x)[[1]]
  if(any(sapply(internals, function(x) length(x$internals))<1))
    stop("plot methods requires internals, set keep=TRUE on twdtwApply() call")
  ts = getTimeSeries(x)[[1]]
  patterns = getPatterns(x)
  
  # Get cost matrix
  df.m = do.call("rbind", lapply(y, function(p){
    tx = index(ts)
    ty = index(shiftDates(patterns[[p]], year=2005))
    m = internals[[p]]$internals[[matrix.name]]
    res = melt(m)
    names(res) = c("Var1","Var2","value")
    res$Pattern = p
    res$tx = tx
    res$ty = ty
    res
  }))
  
  ## Set axis breaks and labels 
  x.labels = pretty_breaks()(range(df.m$tx, na.rm = TRUE))
  timeline = unique( c(df.m$tx, x.labels) )
  x.breaks = zoo( c(unique(df.m$Var2), rep(NA, length(x.labels))), timeline )
  x.breaks = na.approx(x.breaks, rule = 2)
  x.axis = data.frame(x.breaks=x.breaks[x.labels], x.labels = names(x.labels))
  
  fact = 0 
  for(i in seq_along(y)[-1]) fact[i] = fact[i-1] + max(df.m$Var1[df.m$Pattern==y[i]])
  df.m$Var3 = df.m$Var1 + unlist(lapply(seq_along(y), function(i) rep(fact[i], length(which(df.m$Pattern==y[i])) )))
  
  y.axis = do.call("rbind", lapply(y, function(p){
    df = df.m[df.m$Pattern==p,]
    y.labels = pretty_breaks()(range(df$ty, na.rm = TRUE))
    timeline = unique( c(df$ty, y.labels) )
    y.breaks = zoo( c(unique(df$Var3), rep(NA, length(y.labels))), timeline )
    y.breaks = na.approx(y.breaks, rule = 2)
    y.breaks = y.breaks[y.labels]
    data.frame(y.breaks, y.labels=names(y.labels))
  }))
  
  # Plot
  gp = ggplot(data=df.m, aes_string(y='Var3', x='Var2')) +
    facet_wrap(~Pattern, scales = "free", ncol=1) + 
    geom_raster(aes_string(fill='value')) + 
    scale_fill_gradientn(name = legend_name, colours = gray.colors(100, start = 0.1, end = 1)) +
    scale_x_continuous(expand = c(0, 0), breaks=x.axis$x.breaks, labels=x.axis$x.labels) +
    scale_y_continuous(expand = c(0, 0), breaks=y.axis$y.breaks, labels=y.axis$y.labels) +
    xlab("Time series") + 
    ylab("Pattern")
  
  gp
  
}

