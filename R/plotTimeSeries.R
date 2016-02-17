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


#' @title Plotting temporal patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the temporal patterns.
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object, a \code{\link[zoo]{zoo}},
#' or a list of \code{\link[zoo]{zoo}} objects.
#' 
#' @param y An integer. If x is a list \code{y} is an integer or character 
#' (list name(s)). If missing, all elements in the list will be plotted 
#' (up to a maximum of 6).
#' 
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotMatches}}, 
#' \code{\link[dtwSat]{plotClassification}}, and 
#' \code{\link[dtwSat]{plotPatterns}}.
#'  
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun, 
#'         normalize.patterns=TRUE, patterns.length=23)
#'        
#' gp1 = plotTimeSeries(matches)
#' gp1
#' 
#' gp2 = plotTimeSeries(patterns.list)
#' gp2
#' 
#' gp3 = plotTimeSeries(example_ts.list)
#' gp3
#' 
#' @export
plotTimeSeries = function(x, y){
  
  if(!is(x, "list")) x = list(x)
  
  if(missing(y))
    y = names(x)
  
  if(is.null(y)){
    y = paste("Time series",seq_along(x))
    names(x) = y
  } 
  
  if(length(y)>6){
    y = y[1:6]
    x = x[1:6]
  }
    
  x = lapply(x, function(x){
      if(is(x, "twdtw")) x = getTimeSeries(x)
      x
  })
  
  if(any(!sapply(x[y], is.zoo)))
    stop("x is not zoo or list of zoo objects")
  
  # Build data.frame
  df.p = do.call("rbind", lapply(y, function(p)
    data.frame(Time=index(x[[p]]), x[[p]], Series=p)
  ))
  df.p = melt(df.p, id.vars=c("Time","Series"))
  
  # Plot temporal patterns
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    facet_wrap(~Series, scales = "free_x") 
  
  gp
  
}

