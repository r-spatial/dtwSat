#' dtwSat class
#'
#' Class for Multidimensional Time-Weighted DTW results 
#'
#'
#'@section Slots :
#' \describe{
#' \item{\code{call}:}{Object of class \code{call} see 
#' \code{\link[base]{match.call}}}
#' \item{\code{alignments}:}{Object of class \code{list}. Each node of 
#' the list has one attribute, whose lengths are identical to the number of 
#' alignments. \code{from} the alignment starts, \code{to} the alignment ends, 
#' \code{distance} the DTW distance, \code{normalizedDistance} the 
#' normalized DTW distance}
#' \item{\code{mapping}:}{Object of class \code{list}. Each node of the list has 
#' the matching points of the query to the template time series.}
#' \item{\code{internals}:}{Object of class \code{list} similar to 
#' \code{\link[dtw]{dtw}} class from package \pkg{dtw}}
#' }
#'
#' @name dtwSat
#' @aliases dtwSat-class
#' @exportClass dtwSat
#' @author Victor Maus
dtwSat = setClass(
  Class = "dtwSat",
  slots = c(
    call = "call",
    alignments = "list",
    mapping = "list",
    internals = "list"
  ),
  validity = function(object){
    if(!is(object@call, "call")){
      stop("[dtwSat: validation] Invalid call object, class different from call.")
    }else{}
    if(!is(object@alignments, "list")){
      stop("[dtwSat: validation] Invalid TWDTW alignments, class different from data.frema.")
    }else{}
    if(!is(object@mapping, "list")){
      stop("[dtwSat: validation] Invalid TWDTW mapping, class different from list.")
    }else{}
    if(!is(object@internals, "list")){
      stop("[dtwSat: validation] Invalid DTW internals, class different from dtw.")
    }else{}
    return(TRUE)
  }
)


setMethod("initialize",
  signature = "dtwSat",
  definition = 
    function(.Object, call, internals, alignments, mapping){
      .Object@call = new("call")
      .Object@alignments = list(
                                from=numeric(0), 
                                to=numeric(0), 
                                distance=numeric(0), 
                                normalizedDistance=numeric(0))
      .Object@mapping =    list(
                                index1 = numeric(0), 
                                index2 = numeric(0))
      .Object@internals =  list()
      if(!missing(call))
        .Object@call = call
      if(!missing(alignments))
        .Object@alignments = alignments
      if(!missing(mapping))
        .Object@mapping = mapping
      if(!missing(internals))
        .Object@internals = internals
      validObject(.Object)
      return(.Object)
  }
)

setMethod("show", 
          signature = signature(object="dtwSat"),
          definition = function(object){
            cat("Time-Weighted DTW alignment object\n")
            cat("Alignments:\n")
            print(data.frame(object@alignments))
            invisible(NULL)
          }
)

#' @title Plotting Time-Weighted DTW 
#' 
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{plotPath}} and \code{\link[dtwSat]{plotAlignment}}
#' @param x \code{\link[dtwSat]{twdtw}} object
#' @param type A character for the plot type, "path" or "alignment". 
#' Default is "path"
#' @importFrom graphics plot
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
                   image(t(x@internals$costMatrix), col=terrain.colors(100))
            )
          } 
)

#' @title Plotting Time-Weighted DTW 
#' 
#' @description Methods for plotting Time-Weighted DTW objects 
#' returned by twdtw.
#' 
#' @param x \code{\link[dtwSat]{twdtw}} object
#' @param xlab A character, x axis label. Default is "Pattern"
#' @param ylab A character, y axis label. Default is "Time series"
#' @param legend.title A character, the legend for the colour map. 
#' Default is "Warp cost"
#' @param map.colours Colours for the matrix values. 
#' Default is terrain.colors(100)
#' @param x.breaks see y.breaks
#' @param y.breaks control the breaks in the guide, see 
#' \code{\link[ggplot2]{continuous_scale}} in package \pkg{ggplot2}
#' @param x.labels see y.labels
#' @param y.labels a character vector the same length as breaks, see 
#' \code{\link[ggplot2]{continuous_scale}} in package \pkg{ggplot2}
#' @param normalize Normalized distance. Default is FALSE
#' @param show.dist Display dtw distance for each alignment 
#' @param ... additional arguments passed to plotting functions 
#' @docType methods
#' @examples
#' #alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)
#' #plot(alig)
#' @export
plotPath = function(x, ylab="Pattern", xlab="Time series", 
                    legend.title="Warp cost", map.colours=terrain.colors(100),
                    x.breaks=NULL, y.breaks=NULL, x.labels=NULL, y.labels=NULL, 
                    normalize=FALSE, show.dist=FALSE){
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
  # if(is.null(x.breaks)|is.null(x.labels)){
  #   r_x = diff(range(df$Var2)) / as.numeric(diff(range(tx[df$Var2])))
  #   x.labels = pretty.breaks()(range(tx[df$Var2]))
  #   x.breaks = r_x*abs(as.numeric(tx[df$Var2][1] - x.labels))
  # }
  # if(is.null(y.breaks)|is.null(y.labels)){
  #   r_y = diff(range(df$Var1)) / as.numeric(diff(range(ty[df$Var1])))
  #   y.labels = pretty.breaks()(range(ty[df$Var1]))
  #   y.breaks = r_y*abs(as.numeric(ty[df$Var1][1] - y.labels))
  # }
  if(is.null(x.breaks)|is.null(x.labels)){
    x.labels = pretty_breaks()(range(tx[df.m$X2]))
    timeline = unique( c(tx[df.m$X2], x.labels) )
    x.breaks = zoo( c(unique(df.m$X2), rep(NA, length(x.labels))), timeline )
    x.breaks = na.approx(x.breaks)
    x.breaks = x.breaks[x.labels]
  }
  if(is.null(y.breaks)|is.null(y.labels)){
    y.labels = pretty_breaks()(range(ty[df.m$X1]))
    timeline = unique( c(ty[df.m$X1], y.labels) )
    y.breaks = zoo( c(unique(df.m$X1), rep(NA, length(y.labels))), timeline )
    y.breaks = na.approx(y.breaks)
    y.breaks = y.breaks[y.labels]
  }
  
  ## Plot matrix and paths
  gp = ggplot(data=df.m, aes(y=X1, x=X2)) +
    geom_raster(aes(fill=value)) + 
    scale_fill_gradientn(name =legend.title, colours = map.colours) +
    geom_path(data=df.path, aes(y=index1, x=index2, group=alignment)) + 
    scale_y_continuous(expand = c(0, 0), breaks=y.breaks, labels=names(y.labels)) +
    scale_x_continuous(expand = c(0, 0), breaks=x.breaks, labels=names(x.labels)) +
    xlab(xlab) + 
    ylab(ylab)
  
  # Show distance
  if(show.dist){
    text.label = format(x@alignments$distance, digits=2, nsmall = 2)
    if(normalize)
      text.label = format(x@alignments$normalizedDistance, digits=2, nsmall = 2)
    text.x = x@alignments$to + 2
    text.y = rep(max(df.path$index1), length(x))
    gp = gp + annotate("text", x=text.x, y=text.y, label=text.label)
  }
  gp
}






#' @title Plotting Time-Weighted DTW 
#' 
#' @description Methods for plotting Time-Weighted DTW objects 
#' returned by twdtw.
#' 
#' @param x \code{\link[dtwSat]{twdtw}} object
#' @param xlab A character, x axis label
#' @param ylab A character, y axis label
#' @param show.dist Display dtw distance for each alignment 
#' @param ... additional arguments passed to plotting functions 
#' @docType methods
#' @examples
#' #alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50)
#' #plot(alig)
#' @export
plotAlignment = function(x, ylab="Pattern", xlab="Time series", show.dist=FALSE, ...){
  if(length(x@internals)==0)
    stop("plot.dtwSat requires dtwSat internals (set keep=TRUE on twdtw() call)")
  
  m = t(x@internals$costMatrix)
  tx = index(x@internals$template)
  ty = index(x@internals$query)
  
  image(m, col=terrain.colors(100), x=tx, y=ty, xlab=xlab, ylab=ylab, ...)
  contour(m, x=tx, y=ty, add=TRUE)
  for(i in 1:length(x@mapping)){
    lines(tx[x@mapping[[i]]$index2], ty[x@mapping[[i]]$index1], col="red", lwd=2)
    #               if(show.dist)
    #                 text(tail(tx[x$mapping[[i]]$index2],1) , tail(ty[x$mapping[[i]]$index1], 1), round(x$alignments$distance[[i]],2))
  }
}
