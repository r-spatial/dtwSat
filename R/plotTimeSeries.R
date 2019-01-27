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


#' @title Plotting time series 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the temporal patterns.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwTimeSeries}}, 
#' \code{\link[zoo]{zoo}}, or list of class \code{\link[zoo]{zoo}}.
#' @param labels A vector with labels of the time series. If missing, all 
#' elements in the list will be plotted (up to a maximum of 16).
#' @param attr An \link[base]{integer} vector or \link[base]{character} vector 
#' indicating the attribute for plotting. If not declared the function will plot 
#' all attributes.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}} and 
#' \code{\link[dtwSat]{plotPatterns}}
#'  
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' plotTimeSeries(ts)
#' plotTimeSeries(ts, attr="evi")
#' 
#' @export
plotTimeSeries = function(x, labels=NULL, attr){
  
  if(is(x, "twdtwMatches")) x = subset(x@timeseries, labels)
  if(is(x, "twdtwTimeSeries")) x = subset(x, labels) 
  if(is.null(labels)) labels = labels(x)
  new_labels = labels(x)
  labels_tabel = table(new_labels)
  if(any(labels_tabel>1))
    for(p in names(labels_tabel)){
      i = p==labels(x)
      new_labels[i] = paste(new_labels[i], 1:labels_tabel[p])
    }
  x = twdtwTimeSeries(x@timeseries, labels=new_labels)
  labels = new_labels
  if(length(labels)>16) labels = labels[1:16]
      
  # Build data.frame
  if(missing(attr)) attr = names(x[[1]])
  df.p = do.call("rbind", lapply(as.list(x), function(xx){
    ts = xx[[1]][,attr,drop=FALSE]
    data.frame(Time=index(ts), ts, Series=labels(xx)[1])
  }))
  df.p = melt(df.p, id.vars=c("Time","Series"))
  
  # Plot time series 
  gp = ggplot(df.p, aes_string(x="Time", y="value", colour="variable") ) + 
    geom_line() + 
    theme(legend.position = "bottom") + 
    facet_wrap(~Series, scales = "free_x", ncol=1) + 
    guides(colour = guide_legend(title = "Bands")) + 
    ylab("Value") 
  
  gp
  
}


