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
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object
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
#' gp1 = plotPath(alig, n.alignments=1:4)
#' 
#' gp1
#' 
#' gp2 = plotPath(alig, p.names=c("Cotton","Maize"), n.alignments=1:4)
#' 
#' gp2
#' 
#' @export
plotPath = function(x, p.names, n.alignments=NULL){
  
  if(missing(p.names)) {
    p.names = getPatternNames(x)
  } else {
    p.names = getPatternNames(x, p.names)
  }
  
  # Get cost matrix
  df.m = do.call("rbind", lapply(p.names, function(p){
    
    ## Get data
    internals  = getInternals(x, p)
    if(is.null(internals))
      stop("plot methods requires twdtw internals, set keep=TRUE on twdtw() call")
    matching   = getMatches(x, p)
    
    tx = index(internals[[p]]$x)
    ty = index(shiftDate(x=internals[[p]]$pattern, year=2005))
    m = internals[[p]]$costMatrix
    res = melt(m)
    res$Pattern = p
    res$tx = tx
    res$ty = ty
    res
  }))
  
  # Get minimun cost paths
  df.path = do.call("rbind", lapply(p.names, function(p){
    
    ## Get data
    matching = getMatches(x, p)
    
    if(is.null(n.alignments)) n.alignments = seq_along(matching[[p]])
    k = which(n.alignments > length(matching[[p]]))
    if(length(k)>0){
      warning("parameter n.alignments out of bounds")
      n.alignments = n.alignments[-k]
    }
    
    res = do.call("rbind", lapply(n.alignments, function(i){
      data.frame(matching[[p]][[i]], alignment=i)
    }))
    
    res$Pattern = p 
    res 
  }))
  
  ## Set axis breaks and labels 
  x.labels = pretty_breaks()(range(df.m$tx, na.rm = TRUE))
  timeline = unique( c(df.m$tx, x.labels) )
  x.breaks = zoo( c(unique(df.m$Var2), rep(NA, length(x.labels))), timeline )
  x.breaks = na.approx(x.breaks)
  x.axis = data.frame(x.breaks=x.breaks[x.labels], x.labels = names(x.labels))
  
  fact = 0 
  for(i in seq_along(p.names)[-1]) fact[i] = fact[i-1] + max(df.m$Var1[df.m$Pattern==p.names[i]])
  df.m$Var3 = df.m$Var1 + unlist(lapply(seq_along(p.names), function(i) rep(fact[i], length(which(df.m$Pattern==p.names[i])) )))
  df.path$Var3 = df.path$index1 + unlist(lapply(seq_along(p.names), function(i) rep(fact[i], length(which(df.path$Pattern==p.names[i])) )))
  
  y.axis = do.call("rbind", lapply(p.names, function(p){
    df = df.m[df.m$Pattern==p,]
    y.labels = pretty_breaks()(range(df$ty, na.rm = TRUE))
    timeline = unique( c(df$ty, y.labels) )
    y.breaks = zoo( c(unique(df$Var3), rep(NA, length(y.labels))), timeline )
    y.breaks = na.approx(y.breaks)
    y.breaks = y.breaks[y.labels]
    data.frame(y.breaks, y.labels=names(y.labels))
  }))
  
  # Plot
  gp = ggplot(data=df.m, aes_string(y='Var3', x='Var2')) +
    facet_wrap(~Pattern, scales = "free", ncol=1) + 
    geom_raster(aes_string(fill='value')) + 
    scale_fill_gradientn(name = 'Warp cost', colours = terrain.colors(100)) + 
    geom_path(data=df.path, aes_string(y='Var3', x='index2', group='alignment')) + 
    scale_x_continuous(expand = c(0, 0), breaks=x.axis$x.breaks, labels=x.axis$x.labels) +
    scale_y_continuous(expand = c(0, 0), breaks=y.axis$y.breaks, labels=y.axis$y.labels) +
    xlab("Time series") + 
    ylab("Pattern")
  
  gp
  
}

