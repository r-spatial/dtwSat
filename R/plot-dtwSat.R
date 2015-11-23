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
#### twdtw PLOT METHODS


#' @title Plotting twdtw objects
#' 
#' @description Methods for plotting the results of the 
#' Time-Weighted DTW analysis
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object
#' @param type A character for the plot type, ''path'', ''match'', 
#' ''alignment'', ''group'', ''cost''. Default is ''path''
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{twdtw-class}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}}, \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @return A \link[gtable]{gtable} object for methods: ''path'' and 
#' ''cost'', or a \link[ggplot2]{ggplot} object for methods: ''match'', 
#' ''alignment'', and ''group''
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @examples 
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' # Plot paths
#' gp1 = plot(alig, type="path", n.alignments=1:4, show.dist = TRUE)
#' grid.arrange(gp1)
#' 
#' # Plot matches 
#' gp2 = plot(alig, type="match", attr=1)
#' gp2
#' 
#' # Plot alignments 
#' gp3 = plot(alig, type="alignment", attr=c("ndvi","evi"), threshold=4)
#' gp3
#' 
#' ## Plot classification
#' gp4 = plot(alig, type="group", attr="evi", from=as.Date("2009-09-01"),  
#'            to=as.Date("2013-09-01"), by = "6 month", overlap=.3)
#' gp4
#' 
#' # Plot cost matrix
#' gp5 = plot(alig, type="cost", matrix.name="costMatrix")
#' grid.arrange(gp5)
#' 
#' @export
setMethod("plot", 
          signature(x = "twdtw"),
          function(x, type="path", ...){
            if(!is(x,"twdtw"))
              stop("x is not a twdtw object.")
            if(length(getInternals(x))==0)
              stop("plot method requires twdtw internals (set keep.internals=TRUE on dtw() call)")
            pt = pmatch(type,c("path","match","alignment","group","cost"))
            switch(pt,
                   plotPath(x, ...),
                   plotMatch(x, ...),
                   plotAlignment(x, ...),
                   plotGroup(x, ...),
                   plotCostMatrix(x, ...)
            )
          } 
)


#' @title Plotting paths of twdtw alignments
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the minimum paths in the 
#' cost matrix of Time-Weighted DTW
#' 
#' @param x An \code{\link[dtwSat]{twdtw-class}} object
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' will plot the paths for all patterns 
#' @param n.alignments An \link[base]{integer} vector. The alignment indices 
#' to plot. If not declared the function will plot all possible alignments 
#' @param show.dist Display dtw distance for each alignment 
#' @param shift A vector of length 2. These values shift the position
#' of the text in \code{x} and \code{y}, respectively. Argument used 
#' with show.dist
#' @docType methods
#' 
#' @return A \link[gtable]{gtable} object
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#'  
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#'        
#' gp = plotPath(alig, n.alignments=1:4)
#' 
#' grid.arrange(gp)
#' 
#' @export
plotPath = function(x, p.names, n.alignments=NULL, show.dist=FALSE, shift=c(-4,-1)){
  
  if(missing(p.names)) {
    p.names = getPatternNames(x)
  } else {
    p.names = getPatternNames(x, p.names)
  }
  
  last.pattern = tail(p.names, 1)
    
  gp.list = lapply(p.names, function(p){

    ## Get data
    internals  = getInternals(x, p)
    if(is.null(internals))
      stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
    matching   = getMatches(x, p)

    tx = index(internals[[p]]$x)
    ty = index(internals[[p]]$pattern)
    m = internals[[p]]$costMatrix
    df.m = melt(m)
    
    if(is.null(n.alignments)) n.alignments = seq_along(matching[[p]])
    k = which(n.alignments > length(matching[[p]]))
    if(length(k)>0){
      warning("parameter n.alignments out of bounds")
      n.alignments = n.alignments[-k]
    }
    
    df.path = do.call("rbind", lapply(n.alignments, function(i){
      data.frame(matching[[p]][[i]], alignment=i)
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
      ggtitle(p) + 
      xlab("") + 
      ylab("Pattern")
    
    if(p==last.pattern)
      gp = gp + xlab("Time series")
    
    # Show distance
    if(show.dist){
      alig = getAlignments(x, p)
      text.label = format(alig$distance[n.alignments], digits=2, nsmall = 2)
      text.x = unlist(lapply(n.alignments, function(i) tail(matching[[p]][[i]]$index2,1) + shift[1] ) )
      text.y = rep(max(df.path$index1, na.rm = TRUE), length(x)) + shift[2]
      gp = gp + annotate("text", x=text.x, y=text.y, label=text.label)
    }
    gp
  })

  if(length(gp.list)==1)
    return(gp.list[[1]])
  # grid.arrange(grobs = gp.list, nrow=length(p.names))
  arrangeGrob(grobs = gp.list, nrow=length(p.names))
}



#' @title Plotting matching points of twdtw alignments
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the matching points of 
#' Time-Weighted DTW alignments
#' 
#' @param x An \code{\link[dtwSat]{twdtw-class}} object
#' @param p.names A \link[base]{character} or \link[base]{integer}
#' vector with the patterns identification. If not declared the function 
#' will plot one alignment for each pattern in the 
#' \code{\link[dtwSat]{twdtw-class}} object
#' @param n An \link[base]{integer} vector the same length as \code{p.names}.
#' The indices of the alignments for plotting. The alignments in the 
#' \code{\link[dtwSat]{twdtw-class}} object are ordered by TWDTW distance in 
#' ascending order. If not declared the function will plot the best alignment 
#' for each pattern
#' @param attr An \link[base]{integer} or \link[base]{character} vector 
#' indicating the attribute for plotting, \emph{i.e.} a column of the \code{pattern}. 
#' Default is 1
#' @param shift A number, it shifts the pattern position in the \code{x}
#' direction. Default is 0.5
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' gp = plotMatch(alig)
#' gp
#' 
#' gp = plotMatch(x=alig, p.names=1)
#' gp
#' 
#' gp = plotMatch(x=alig, p.names="Cotton")
#' gp
#' 
#' gp = plotMatch(x=alig, p.names=c("Cotton","Cotton","Cotton","Cotton"), 
#'      n = c(1:4))
#' gp
#' 
#' @export
plotMatch = function(x, p.names, n, attr=1, shift=0.5){
  
  if(missing(p.names)) {
    p.names = getPatternNames(x)
  } else {
    p.names = getPatternNames(x, p.names)
  }
  
  if(missing(n)) n = rep(1, length(p.names))
  if(!length(n)==length(p.names))
    stop("n is not the same length as p.names")
  
  names(n) = p.names
 
  ## Get data
  internals  = getInternals(x, p.names)
  if(is.null(internals))
    stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
  matching   = getMatches(x, p.names)
  alignments = getAlignments(x, p.names)
  
  xx = internals[[p.names[1]]]$x[,attr,drop=FALSE]
  tx = index(xx)

  y.labels = pretty_breaks()(range(xx, na.rm = TRUE))
  y.breaks = y.labels
  
  # Get time series 
  df.x = data.frame(Time=tx, xx)
  
  # Get matching points
  df.list = lapply(seq_along(p.names), function(i){
    p = p.names[i]
    yy = internals[[p]]$pattern[,attr,drop=FALSE]
    ty = index(yy)
    
    map = data.frame(matching[[p]][[n[i]]])
    delay = tx[map$index2[1]]-ty[1]
    if(delay>0)
      delay = delay + diff(range(ty, na.rm = TRUE))*shift
    if(delay<0)
      delay = delay - diff(range(ty, na.rm = TRUE))*shift
    
    df.pt = data.frame(Time=ty[map$index1]+delay, yy[map$index1,,drop=FALSE]+max(xx, na.rm = TRUE))
    df.match.pt = df.pt
    df.match.pt$alig = paste(1:nrow(map),p,n[i],sep="_")
    df.match.x = df.x[map$index2,]
    df.match.x$alig = paste(1:nrow(map),p,n[i],sep="_")
    df.match = rbind(df.match.pt, df.match.x)
    df.pt$Pattern = paste(p,n[i])
    list(match=df.match, pt=df.pt)
  })
  
  df.pt = do.call("rbind", lapply(df.list, function(df) df$pt))
  df.match = do.call("rbind", lapply(df.list, function(df) df$match))
  
  attr_names = names(df.x)[2]
  gp = ggplot(data=df.x, aes_string(x='Time', y=eval(attr_names))) +
    geom_line() +
    geom_line(data=df.pt, aes_string(x='Time', y=eval(attr_names), 
                                     group='Pattern', colour='Pattern')) + 
    geom_line(data=df.match, linetype = 2, colour = "grey", 
              aes_string(x='Time', y=eval(attr_names), group='alig')) + 
    scale_y_continuous(breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(attr_names) 
  
  gp
  
}



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



#' @title Plotting cost matrices
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the internal matrices
#' 
#' @param x An \code{\link[dtwSat]{twdtw-class}} object
#' @param matrix.name A character. The name of the matrix to plot,
#' "costMatrix", "localMatrix", or "timeWeight". 
#' Default is "costMatrix"
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' will plot the matrices for all patterns
#' 
#' @return A \link[gtable]{gtable} object 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' # Plot local cost for all patterns
#' gp = plotCostMatrix(x=alig, matrix.name="localMatrix")
#' grid.arrange(gp)
#' 
#' # Plot time weight for all patterns
#' gp = plotCostMatrix(x=alig, matrix.name="timeWeight")
#' grid.arrange(gp)
#' 
#' # Plot accumulated cost for all patterns
#' gp = plotCostMatrix(x=alig, matrix.name="costMatrix")
#' grid.arrange(gp)
#' 
#' # Plot accumulated cost for Cotton and Soybean
#' gp = plotCostMatrix(x=alig, matrix.name="costMatrix", 
#'      p.names=c("Cotton","Soybean"))
#' grid.arrange(gp)
#' 
#' 
#' @export
plotCostMatrix = function(x, matrix.name="costMatrix", p.names){
  
  pt = pmatch(matrix.name,c("costMatrix", "localMatrix", "timeWeight"))
  if(is.na(pt))
    stop("matrix.name is not costMatrix, localMatrix, or timeWeight")
  
  if(missing(p.names)) {
    p.names = getPatternNames(x)
  } else {
    p.names = getPatternNames(x, p.names)
  }
  
  last.pattern = tail(p.names, 1)
  
  gp.list = lapply(p.names, function(p){
    
    ## Get data
    internals  = getInternals(x, p)
    if(is.null(internals))
      stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
    
    tx = index(internals[[p]]$x)
    ty = index(internals[[p]]$pattern)
    m = internals[[p]][[matrix.name]]
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
      ggtitle(p) + 
      xlab("") + 
      ylab("Pattern")
    
    if(p==last.pattern)
      gp = gp + xlab("Time series")
    
    gp
    
  })
  
  if(length(gp.list)==1)
    return(gp.list[[1]])
  # grid.arrange(grobs = gp.list, nrow=length(p.names))
  arrangeGrob(grobs = gp.list, nrow=length(p.names))
 
}
