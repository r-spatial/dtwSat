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
#' @param x An \code{\link[dtwSat]{twdtw-class}} object, a \code{\link[zoo]{zoo}},
#' or a list of \code{\link[zoo]{zoo}} objects.
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' will plot the paths for all patterns. 
#' @param year An integer. The base year to shift the dates of the time series to. 
#' If NULL then it does not shif the time series. Default is 2005. 
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
#' \code{\link[dtwSat]{plotMatches}}, and 
#' \code{\link[dtwSat]{plotClassification}}.
#'  
#' @examples
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun, 
#'         normalize.patterns=TRUE, patterns.length=23)
#'        
#' gp1 = plotPatterns(matches)
#' gp1
#' 
#' gp2 = plotPatterns(patterns.list)
#' gp2
#' 
#' gp3 = plotPatterns(patterns.list$Soybean)
#' gp3
#' 
#' @export
plotPatterns = function(x, p.names, year=2005){
  
  # Get temporal patterns
  if(is(x, "twdtw")){
    if(missing(p.names)) {
      p.names = getPatternNames(x)
    } else {
      p.names = getPatternNames(x, p.names)
    }
    x = getPatterns(x, p.names)
  }
  
  if(is(x, "zoo")) x = list(x)
  
  if(missing(p.names))
    p.names = names(x)
  
  if(is.null(p.names)){
    p.names = paste("Pattern",seq_along(x))
    names(x) = p.names
  } 
  
  if(any(!sapply(x[p.names], is.zoo)))
    stop("patterns should be a list of zoo objects")
  
  # Shift dates 
  if(!is.null(year)) x = lapply(x[p.names], shiftDates, year=year)
  
  # Build data.frame
  df.p = do.call("rbind", lapply(p.names, function(p)
    data.frame(Time=index(x[[p]]), x[[p]], Pattern=p)
  ))
  df.p = melt(df.p, id.vars=c("Time","Pattern"))
  
  # Plot temporal patterns
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    facet_wrap(~Pattern) + 
    scale_x_date(labels = date_format("%b"))
  
  gp
  
}

