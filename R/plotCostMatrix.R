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
#' @return A \link[ggplot2]{ggplot} object 
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
#' gp1 = plotCostMatrix(x=alig, matrix.name="localMatrix")
#' gp1
#' 
#' # Plot time weight for all patterns
#' gp2 = plotCostMatrix(x=alig, matrix.name="timeWeight")
#' gp2
#' 
#' # Plot accumulated cost for all patterns
#' gp3 = plotCostMatrix(x=alig, matrix.name="costMatrix")
#' gp3
#' 
#' # Plot accumulated cost for Cotton and Soybean
#' gp4 = plotCostMatrix(x=alig, matrix.name="costMatrix", 
#'      p.names=c("Cotton","Soybean"))
#' 
#' gp4
#' 
#' @export
plotCostMatrix = function(x, matrix.name="costMatrix", p.names){
  
  pt = pmatch(matrix.name,c("costMatrix", "localMatrix", "timeWeight"))
  if(is.na(pt))
    stop("matrix.name is not costMatrix, localMatrix, or timeWeight")
  
  legend_name = c("Warp cost", "Local cost", "Time weight")[pt]
  
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
    m = internals[[p]][[matrix.name]]
    res = melt(m)
    res$Pattern = p
    res$tx = tx
    res$ty = ty
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
    scale_fill_gradientn(name = legend_name, colours = gray.colors(100, start = 0.1, end = 1)) +
    scale_x_continuous(expand = c(0, 0), breaks=x.axis$x.breaks, labels=x.axis$x.labels) +
    scale_y_continuous(expand = c(0, 0), breaks=y.axis$y.breaks, labels=y.axis$y.labels) +
    xlab("Time series") + 
    ylab("Pattern")
  
  gp
  
}

