###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Münster (WWU), Germany                  #
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
#' ''classify''
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}}, and 
#' \code{\link[dtwSat]{plotClassify}}
#' Default is "path"
#' 
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso \code{\link[dtwSat]{dtwSat-class}}, \code{\link[dtwSat]{dtwSat}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}}, and 
#' \code{\link[dtwSat]{plotClassify}}
#' 
#' @examples 
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, alignments = 4, keep = TRUE)
#' 
#' ## Plot path
#' gp1 = plot(alig, type="path", normalize=TRUE, show.dist = TRUE)
#' gp1
#' 
#' ## Plot alignment
#' gp2 = plot(alig, type="alignment", attribute="evi", alignment=1, shift=0.5)
#' gp2
#' 
#' ## Plot classification
#' gp = plot(alig, type="classify", attribute="evi", 
#'          from=as.Date("2009-09-01"),  to=as.Date("2013-09-01"), 
#'          by = "6 month", normalized=TRUE, overlap=.7)
#' gp
#' 
#' @export
setMethod("plot", 
          signature(x = "dtwSat"),
          function(x, type="path", ...){
            if(!is(x,"dtwSat"))
              stop("x is not a dtwSat object.")
            if(length(getInternals(x))==0)
              stop("plot method requires dtwSat internals (set keep.internals=TRUE on dtw() call)")
            pt = pmatch(type,c("path","alignment","classify"))
            switch(pt,
                   plotPath(x, ...),
                   plotAlignment(x, ...),
                   plotClassify(x, ...)
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
#' @param normalize Plot normalized distance. Default is FALSE
#' @param alignments A vector. The alignment indices to 
#' plot. NULL will plot all possible alignments 
#' @param show.dist Display dtw distance for each alignment 
#' @param shift A vector of length 2. These values shift the position
#' of the text in \code{x} and \code{y}, respectively. Argument used 
#' with show.dist.
#' @docType methods
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' #' @seealso  \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{plotAlignment}}, 
#' and \code{\link[dtwSat]{plotClassify}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, alignments = 4, keep = TRUE)
#' gp = plotPath(alig, normalize=TRUE, show.dist = TRUE)
#' gp
#' 
#' @export
plotPath = function(x, normalize=FALSE, alignments=NULL, 
                    show.dist=FALSE, shift=c(-4,-1)){
  ## Get data
  internals = getInternals(x)
  mapping = getMatches(x)
  alig = getAlignments(x)
  
  tx = index(internals$template)
  ty = index(internals$query)
  m = internals$costMatrix
  df.m = melt(m)
  
  if(normalize)
    df.m$value = df.m$value/nrow(m)
  
  if(is.null(alignments))
    alignments = seq_along(mapping)
  if(any(alignments > length(mapping)))
    stop("parameter alignments out of bounds")

  df.path = do.call("rbind", lapply(alignments, function(i){
    data.frame(mapping[[i]], alignment=i)
  }))

  ## Set axis breaks and labels 
  x.labels = pretty_breaks()(range(tx[df.m$Var2]))
  timeline = unique( c(tx[df.m$Var2], x.labels) )
  x.breaks = zoo( c(unique(df.m$Var2), rep(NA, length(x.labels))), timeline )
  x.breaks = na.approx(x.breaks)
  x.breaks = x.breaks[x.labels]
  y.labels = pretty_breaks()(range(ty[df.m$Var1]))
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
    text.label = format(alig$distance[alignments], digits=2, nsmall = 2)
    if(normalize)
      text.label = format(alig$normalizedDistance[alignments], digits=2, nsmall = 2)
    text.x = unlist(lapply(alignments, function(i) tail(mapping[[i]]$index2,1) + shift[1] ) )
    text.y = rep(max(df.path$index1), length(x)) + shift[2]
    gp = gp + annotate("text", x=text.x, y=text.y, label=text.label)
  }
  gp
}



#' @title Plotting alignments of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the alignments of 
#' Time-Weighted DTW analysis
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param alignment An integer, the index of the alignment to plot
#' @param attribute An integer or a character indicating the attribute 
#' for plotting, \emph{i.e.} the column of the \code{query}. Default is 1
#' @param shift A number, it shifts the pattern position in the \code{x}
#' direction. Default is 0.5
#' @docType methods
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, 
#' \code{\link[dtwSat]{plotPath}}, and 
#' \code{\link[dtwSat]{plotClassify}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, alignments = 4, keep = TRUE)
#' gp = plotAlignment(alig, attribute="evi", alignment=1, shift=0.5)
#' gp
#' 
#' @export
plotAlignment = function(x, alignment, attribute=1, shift=0.5){
  
  ## Get data
  internals = getInternals(x)
  mapping = getMatches(x)
  alignments = getAlignments(x)
  
  tx = index(internals$template)
  ty = index(internals$query)
  
  xx = internals$template[,attribute]
  yy = internals$query[,attribute]
  
  map = data.frame(mapping[[alignment]])
  
  delay = tx[map$index2[1]]-ty[1]
  if(delay>0)
    delay = delay + diff(range(ty))*shift
  if(delay<0)
    delay = delay - diff(range(ty))*shift
  
  y.labels = pretty_breaks()(range(xx))
  y.breaks = y.labels
  
  df.ts = data.frame(Time=tx, value=xx)
  
  df.pt = data.frame(Time=ty[map$index1]+delay, value=yy[map$index1]+max(xx))
  
  df.match.pt = df.pt
  df.match.pt$p = p=1:nrow(map)
  df.match.ts = df.ts[map$index2,]
  df.match.ts$p = p=1:nrow(map)
  df.match = rbind(df.match.pt, df.match.ts)
  if(!is(attribute,"character"))
    attribute="Value"
  gp = ggplot(data=df.ts, aes_string(x='Time', y='value')) +
    geom_line() +
    geom_line(data=df.pt, aes_string(x='Time', y='value'), col="red") + 
    geom_line(data=df.match, linetype = 2, colour = "grey", aes_string(x='Time', y='value', group='p')) + 
    scale_y_continuous(breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(toupper(attribute)) + 
    xlab("Time")
  gp
}



#' @title Plotting classification of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the classification of 
#' intervals
#' 
#' @param x An \code{\link[dtwSat]{dtwSat-class}} object
#' @param attribute An integer or a character indicating the attribute 
#' for plotting, \emph{i.e.} the column of the \code{query}. Default is 1
#' @param ... additional arguments passed to \code{\link[dtwSat]{classfyIntervals}}
#' 
#' @return object of class \link[ggplot2]{ggplot}
#' 
#' @seealso  \code{\link[dtwSat]{dtwSat}}, 
#' \code{\link[dtwSat]{plotPath}}, and \code{\link[dtwSat]{plotAlignment}}
#' 
#' @examples
#' malig = mtwdtw(query.list, template, weight = "logistic", 
#'                alpha = 0.1, beta = 100)
#' 
#' gp = plotClassify(x=malig, attribute="evi", from=as.Date("2009-09-01"),  
#'              to=as.Date("2013-09-01"), by = "6 month",
#'              normalized=TRUE, overlap=.7)
#' gp
#' 
#' @export
plotClassify = function(x, attribute=1, ...){
  
  ## Get data
  internals = getInternals(x)
  best_class = classfyIntervals(x, ...)
  
  tx = index(internals$template)
  
  xx = internals$template[,attribute]
  
  df.ts = data.frame(Time=tx, value=xx)

  y.labels = pretty_breaks()(range(xx))
  y.breaks = y.labels
  df.pol = do.call("rbind", lapply(1:nrow(best_class), function(i){
    data.frame(
      Time = c(best_class$from[i], best_class$to[i], best_class$to[i], best_class$from[i]),
      Group = rep(i, 4),
      Query = rep(best_class$query[i], 4),
      value = rep(range(y.breaks), each=2))
  }))
  df.pol$Group = factor(df.pol$Group)
  df.pol$Query = factor(df.pol$Query)

  if(!is(attribute,"character"))
    attribute="Value"

  gp = ggplot() +
    geom_polygon(data=df.pol, aes_string(x='Time', y='value', 
                                         group='Group', fill='Query'), alpha=.7) +
    scale_fill_brewer(palette="Set3") + 
    geom_line(data=df.ts, aes_string(x='Time', y='value')) +
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(toupper(attribute)) + 
    xlab("Time")
  gp
}

