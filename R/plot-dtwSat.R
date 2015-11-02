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
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### dtwSat PLOT METHODS


#' @title Plotting dtwSat objects
#' 
#' @description Methods for plotting the results of the 
#' Time-Weighted DTW analysis
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @param x A \code{\link[dtwSat]{dtwSat-class}} object
#' @param type A character for the plot type, ''path'', ''alignment'', or 
#' ''classify''. Default is "path"
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}}, and 
#' \code{\link[dtwSat]{plotClassify}}
#' 
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso \code{\link[dtwSat]{dtwSat-class}}, \code{\link[dtwSat]{dtwSat}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}}, and 
#' \code{\link[dtwSat]{plotClassify}}
#' 
#' @examples 
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, n.alignments = 4, keep = TRUE)
#' 
#' ## Plot path
#' gp1 = plot(alig, type="path", show.dist = TRUE)
#' gp1
#' 
#' ## Plot alignment
#' gp2 = plot(alig, type="alignment", n=1, attr="evi", shift=0.5)
#' gp2
#' 
#' ## Plot classification
#' gp3 = plot(alig, type="classify", attr="evi", 
#'          from=as.Date("2009-09-01"),  to=as.Date("2013-09-01"), 
#'          by = "6 month", overlap=.7)
#' gp3
#' 
#' @export
setMethod("plot", 
          signature(x = "dtwSat"),
          function(x, type="path", ...){
            if(!is(x,"dtwSat"))
              stop("x is not a dtwSat object.")
            if(length(getInternals(x))==0)
              stop("plot method requires dtwSat internals (set keep.internals=TRUE on dtw() call)")
            pt = pmatch(type,c("path","alignment","classify","cost"))
            switch(pt,
                   plotPath(x, ...),
                   plotAlignment(x, ...),
                   plotClassify(x, ...),
                   plotCostMatrix(x, ...)
            )
          } 
)


#' @title Plotting paths of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the minimum paths in the 
#' cost matrix of Time-Weighted DTW
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param n.alignments A vector. The alignment indices to 
#' plot. NULL will plot all possible alignments 
#' @param show.dist Display dtw distance for each alignment 
#' @param shift A vector of length 2. These values shift the position
#' of the text in \code{x} and \code{y}, respectively. Argument used 
#' with show.dist
#' @param normalized Use normalized TWDTW distance. Default is TRUE
#' @docType methods
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{plotAlignment}}, 
#' \code{\link[dtwSat]{plotClassify}}, and \code{\link[dtwSat]{plotCostMatrix}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, n.alignments = 4, keep = TRUE)
#' gp = plotPath(alig, show.dist = TRUE)
#' gp
#' 
#' @export
plotPath = function(x, n.alignments=NULL, show.dist=FALSE, shift=c(-4,-1), normalized){
  
  if (!missing(normalized))
    warning("argument normalized is deprecated and is scheduled to be removed in the next version", 
            call. = FALSE)
  
  ## Get data
  internals = getInternals(x)
  mapping = getMatches(x)
  alig = getAlignments(x)
  
  tx = index(internals$timeseries)
  ty = index(internals$query)
  m = internals$costMatrix
  df.m = melt(m)
  
  if(is.null(n.alignments))
    n.alignments = seq_along(mapping)
  if(any(n.alignments > length(mapping)))
    stop("parameter alignments out of bounds")

  df.path = do.call("rbind", lapply(n.alignments, function(i){
    data.frame(mapping[[i]], alignment=i)
  }))

  ## Set axis breaks and labels 
  x.labels = pretty_breaks()(range(tx[df.m$Var2], na.rm = TRUE))
  timeline = unique( c(tx[df.m$Var2], x.labels) )
  x.breaks = zoo( c(unique(df.m$Var2), rep(NA, length(x.labels))), timeline )
  x.breaks = na.approx(x.breaks)
  x.breaks = x.breaks[x.labels]
  y.labels = pretty_breaks()(range(ty[df.m$Var1], na.rm = TRUE))
  timeline = unique( c(ty[df.m$Var1], y.labels) )
  y.breaks = zoo( c(unique(df.m$Var1), rep(NA, length(y.labels))), timeline )
  y.breaks = na.approx(y.breaks)
  y.breaks = y.breaks[y.labels]

  ## Plot matrix and paths
  gp = ggplot(data=df.m, aes_string(y='Var1', x='Var2')) +
    geom_raster(aes_string(fill='value')) + 
    scale_fill_gradientn(name = 'Warp cost', colours = terrain.colors(100)) +
    geom_path(data=df.path, aes_string(y='index1', x='index2', group='alignment')) + 
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=names(y.labels)) +
    scale_x_continuous(expand = c(0, 0), breaks=x.breaks, labels=names(x.labels)) +
    xlab("Time series") + 
    ylab("Pattern")

  # Show distance
  if(show.dist){
    text.label = format(alig$distance[n.alignments], digits=2, nsmall = 2)
    text.x = unlist(lapply(n.alignments, function(i) tail(mapping[[i]]$index2,1) + shift[1] ) )
    text.y = rep(max(df.path$index1, na.rm = TRUE), length(x)) + shift[2]
    gp = gp + annotate("text", x=text.x, y=text.y, label=text.label)
  }
  gp
}



#' @title Plotting alignment of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the alignments of 
#' Time-Weighted DTW analysis
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param n An integer, the index of the alignment to plot
#' @param attr An integer or a character indicating the attribute 
#' for plotting, \emph{i.e.} the column of the \code{query}. Default is 1
#' @param shift A number, it shifts the pattern position in the \code{x}
#' direction. Default is 0.5
#' @docType methods
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{plotPath}},  
#' \code{\link[dtwSat]{plotClassify}}, and \code{\link[dtwSat]{plotCostMatrix}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, n.alignments = 4, keep = TRUE)
#' gp = plotAlignment(alig, n=1, attr="evi", shift=0.5)
#' gp
#' 
#' @export
plotAlignment = function(x, n, attr=1, shift=0.5){
  
  ## Get data
  internals = getInternals(x)
  mapping = getMatches(x)
  alignments = getAlignments(x)
  
  tx = index(internals$timeseries)
  ty = index(internals$query)
  
  xx = internals$timeseries[,attr]
  yy = internals$query[,attr]
  
  map = data.frame(mapping[[n]])
  
  delay = tx[map$index2[1]]-ty[1]
  if(delay>0)
    delay = delay + diff(range(ty, na.rm = TRUE))*shift
  if(delay<0)
    delay = delay - diff(range(ty, na.rm = TRUE))*shift
  
  y.labels = pretty_breaks()(range(xx, na.rm = TRUE))
  y.breaks = y.labels
  
  df.timeseries = data.frame(Time=tx, value=xx)
  
  df.pt = data.frame(Time=ty[map$index1]+delay, value=yy[map$index1]+max(xx, na.rm = TRUE))
  
  df.match.pt = df.pt
  df.match.pt$p = p=1:nrow(map)
  df.match.timeseries = df.timeseries[map$index2,]
  df.match.timeseries$p = p=1:nrow(map)
  df.match = rbind(df.match.pt, df.match.timeseries)
  if(!is(attr,"character"))
    attr="Value"
  gp = ggplot(data=df.timeseries, aes_string(x='Time', y='value')) +
    geom_line() +
    geom_line(data=df.pt, aes_string(x='Time', y='value'), col="red") + 
    geom_line(data=df.match, linetype = 2, colour = "grey", aes_string(x='Time', y='value', group='p')) + 
    scale_y_continuous(breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(toupper(attr)) + 
    xlab("Time")
  gp
}


#' @title Plotting classification of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the classification of 
#' time intervals
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param attr A vector of integer or a character indicating 
#' the attribute for plotting, \emph{i.e.} the column of the \code{query}. 
#' Default is 1
#' @param ... additional arguments passed to \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotAlignment}}, \code{\link[dtwSat]{plotCostMatrix}}, and 
#' \code{\link[dtwSat]{classifyIntervals}}
#' 
#' @examples
#' malig = mtwdtw(query.list, timeseries=template, weight = "logistic", 
#'                alpha = 0.1, beta = 100)
#' 
#' gp = plotClassify(x=malig, attr=c("ndvi","evi"), from=as.Date("2009-09-01"),  
#'              to=as.Date("2013-09-01"), by = "6 month", overlap=.7)
#' gp
#' 
#' @export
plotClassify = function(x, attr=NULL, ...){

  ## Get data
  internals = getInternals(x)
  best_class = classifyIntervals(x, ...)
  tx = index(internals$timeseries)
  
  I = min(best_class$from, na.rm = TRUE)-30 <= index(internals$timeseries) & 
      index(internals$timeseries) <= max(best_class$to, na.rm = TRUE)+30
  xx = internals$timeseries[I,]
  if(!is.null(attr))
    xx = internals$timeseries[I,attr]
  df.timeseries = melt(data.frame(Time=index(xx), xx), id="Time")
  if(is.null(df.timeseries$variable))
    df.timeseries$variable = names(internals$timeseries)[attr]
  if(length(attr)==1 | any(is.na(df.timeseries$variable)))
    df.timeseries$variable = attr
  
  y.labels = pretty_breaks()(range(df.timeseries$value, na.rm = TRUE))
  y.breaks = y.labels
  df.pol = do.call("rbind", lapply(1:nrow(best_class), function(i){
    data.frame(
      Time = c(best_class$from[i], best_class$to[i], best_class$to[i], best_class$from[i]),
      Group = rep(i, 4),
      Query = rep(best_class$query[i], 4),
      value = rep(range(y.breaks, na.rm = TRUE), each=2))
  }))
  df.pol$Group = factor(df.pol$Group)
  df.pol$Query = factor(df.pol$Query)
    
  gp = ggplot() +
    geom_polygon(data=df.pol, aes_string(x='Time', y='value', 
                                         group='Group', fill='Query'), alpha=.7) +
    scale_fill_brewer(palette="Set3") + 
    geom_line(data=df.timeseries, aes_string(x='Time', y='value', colour='variable')) +
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab("Value") + 
    xlab("Time")
  gp
}



#' @title Plotting cost matrices
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the internal matrices
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param matrix.name A character. The name of the matrix to plot,
#' "costMatrix", "localMatrix", or "timeWeight". 
#' Default is "costMatrix"
#' 
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}},
#' and \code{\link[dtwSat]{plotClassify}}
#' 
#' @examples
#' names(query.list)
#' alig = twdtw(query.list[["Soybean"]], timeseries=template, weight = "logistic", 
#'                alpha = 0.1, beta = 100, keep=TRUE)
#' 
#' gp = plotCostMatrix(x=alig, matrix.name="localMatrix")
#' gp
#' 
#' gp = plotCostMatrix(x=alig, matrix.name="timeWeight")
#' gp
#' 
#' gp = plotCostMatrix(x=alig, matrix.name="costMatrix")
#' gp
#' 
#' @export
plotCostMatrix = function(x, matrix.name="costMatrix"){
  
  pt = pmatch(matrix.name,c("costMatrix", "localMatrix", "timeWeight"))
  if(is.na(pt))
    stop("matrix.name is not costMatrix, localMatrix, or timeWeight")
  
  ## Get data
  internals = getInternals(x)

  tx = index(internals$timeseries)
  ty = index(internals$query)
  m = internals[[matrix.name]]
  df.m = melt(m)

  ## Set axis breaks and labels 
  x.labels = pretty_breaks()(range(tx[df.m$Var2], na.rm = TRUE))
  timeline = unique( c(tx[df.m$Var2], x.labels) )
  x.breaks = zoo( c(unique(df.m$Var2), rep(NA, length(x.labels))), timeline )
  x.breaks = na.approx(x.breaks)
  x.breaks = x.breaks[x.labels]
  y.labels = pretty_breaks()(range(ty[df.m$Var1], na.rm = TRUE))
  timeline = unique( c(ty[df.m$Var1], y.labels) )
  y.breaks = zoo( c(unique(df.m$Var1), rep(NA, length(y.labels))), timeline )
  y.breaks = na.approx(y.breaks)
  y.breaks = y.breaks[y.labels]
  
  legend_name = c("Warp cost", "Local cost", "Time weight")[pt]
  
  ## Plot matrix and paths
  gp = ggplot(data=df.m, aes_string(y='Var1', x='Var2')) +
    geom_raster(aes_string(fill='value')) + 
    scale_fill_gradientn(name = legend_name, colours = gray.colors(100, start = 0.1, end = 1)) +
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=names(y.labels)) +
    scale_x_continuous(expand = c(0, 0), breaks=x.breaks, labels=names(x.labels)) +
    xlab("Time series") + 
    ylab("Pattern")
  gp
}
