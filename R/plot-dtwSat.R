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
#' Time-Weighted DTW analsys
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @param x \code{\link[dtwSat]{dtwSat-class}} object
#' @param type A character for the plot type, "path" or "alignment"
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{plotPath}} and \code{\link[dtwSat]{plotAlignment}}
#' Default is "path"
#' 
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}}
#' 
#' @examples 
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 50, alignments = 4, keep = TRUE)
#' 
#' ## Plot path
#' gp1 = plot(alig, type="path", normalize=TRUE, show.dist = TRUE)
#' gp1
#' 
#' ## Plot alignment
#' gp2 = plot(alig, type="alignment", dimension="evi", alignment=1, shift=0.5)
#' gp2
#' 
#' @export
setMethod("plot", 
          signature(x = "dtwSat"),
          function(x, type="path", ...){
            if(!is(x,"dtwSat"))
              stop("x is not a dtwSat object.")
            if(length(x@internals)==0) 
              stop("plot method requires dtwSat internals (set keep.internals=TRUE on dtw() call)")
            pt = pmatch(type,c("path","alignment"))
            switch(pt,
                   plotPath(x, ...),
                   plotAlignment(x, ...)
            )
          } 
)


#' @title Plotting paths of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the paths in the cost matrix 
#' of the \code{\link[dtwSat]{dtwSat-class}} objects.
#' 
#' @param x \code{\link[dtwSat]{dtwSat-class}} object
#' @param normalize Normalized distance. Default is FALSE
#' @param show.dist Display dtw distance for each alignment 
#' @param shift A vector whose values shift the position of the text 
#' in the plot. Argument Used with show.dist.
#' @docType methods
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' #' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotAlignment}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 50, alignments = 4, keep = TRUE)
#' gp = plotPath(alig, normalize=TRUE, show.dist = TRUE)
#' gp
#' 
#' @export
plotPath = function(x, normalize=FALSE, show.dist=FALSE, shift=c(-4,-1)){
  ## Get data
  tx = index(x@internals$template)
  ty = index(x@internals$query)
  m = x@internals$costMatrix
  df.m = melt(m)

  if(normalize)
    df.m$value = df.m$value/nrow(m)
  df.path = do.call("rbind", lapply(seq_along(x@mapping), function(i){
    data.frame(x@mapping[[i]], alignment=i)
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
    text.label = format(x@alignments$distance, digits=2, nsmall = 2)
    if(normalize)
      text.label = format(x@alignments$normalizedDistance, digits=2, nsmall = 2)
    text.x = x@alignments$to + shift[1]
    text.y = rep(max(df.path$index1), length(x)) + shift[2]
    gp = gp + annotate("text", x=text.x, y=text.y, label=text.label)
  }
  gp
}






#' @title Plotting alignments of dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the alignments in the cost 
#' matrix of the \code{\link[dtwSat]{dtwSat-class}} objects.
#' 
#' @param x \code{\link[dtwSat]{dtwSat-class}} object
#' @param alignment An integer, the index of the alignment
#' @param dimension An integer or a character with the name of the 
#' column in query and template to use in the plot
#' @param shift A number, it shifts the pattern position.
#' Default is 0.5
#' @docType methods
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, alignments = 4, keep = TRUE)
#' gp = plotAlignment(alig, dimension="evi", alignment=1, shift=0.5)
#' gp
#' 
#' @export
plotAlignment = function(x, alignment, dimension, shift=0.5){
  
  tx = index(x@internals$template)
  ty = index(x@internals$query)
  
  xx = x@internals$template[,dimension]
  yy = x@internals$query[,dimension]
  
  map = data.frame(x@mapping[[alignment]])
  
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
  
  gp = ggplot(data=df.ts, aes_string(x='Time', y='value')) +
    geom_line() +
    geom_line(data=df.pt, aes_string(x='Time', y='value'), col="red") + 
    geom_line(data=df.match, linetype = 2, colour = "grey", aes_string(x='Time', y='value', group='p')) + 
    scale_y_continuous(breaks=y.breaks, labels=y.labels) +
    scale_x_date(breaks=waiver(), labels=waiver()) +
    ylab(toupper(dimension)) + 
    xlab("Time")
  gp
}




