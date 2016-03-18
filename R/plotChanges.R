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


#' @title Plotting changes 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting changes over time.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwRaster}}.
#' @param time.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the layers to plot. For plot type ''change'' the minimum length 
#' is two.
#' @param time.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the layers. It must have the same 
#' length as time.levels. Default is NULL.
#' @param class.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the levels of the raster values. Default is NULL. 
#' @param class.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the raster values. It must have the same 
#' length as class.levels. Default is NULL.
#' @param class.colors a set of aesthetic values. It must have the same 
#' length as class.levels. Default is NULL. See 
#' \link[ggplot2]{scale_fill_manual} for details.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}},
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[dtwSat]{plotArea}}, 
#' \code{\link[dtwSat]{plotMaps}}, and
#' \code{\link[dtwSat]{plotDistance}}.
#'  
#' @examples
#' \dontrun{
#' 
#' }
#' @export
plotChanges = function(x, time.levels=NULL, time.labels=NULL, class.levels=NULL, class.labels=NULL, class.colors=NULL){
  plot(x, type="changes", time.levels=time.levels, time.labels=time.labels, class.levels=class.levels, class.labels=class.labels, class.colors=class.colors)
}

.plotChanges = function(x, time.levels, time.labels, class.levels, class.labels, class.colors){

  if(length(time.levels)<2)
    stop("the length of time.levels is shorter than two")
  
  df = do.call("rbind", lapply(seq_along(time.levels)[-1], function(l){
    from = raster::subset(x=x, subset=time.levels[l-1])[]
    to   = raster::subset(x=x, subset=time.levels[l]  )[]
    res = data.frame(from, to)
    res = data.frame(xtabs(~ from + to, res) / nrow(res))
    res$layer = paste0(time.labels[l-1],"-",time.labels[l])
    res
  }))
  
  df$from = factor(df$from, levels = class.levels, labels = class.labels)
  df$to   = factor(df$to, levels = class.levels, labels = class.labels)
  df$from = factor(df$from)
  df$to   = factor(df$to)
  I = df$from!=df$to
  
  # Plot change 
  gp = ggplot() +
    geom_bar(data=df[I,], aes_string(x="to", y="Freq", fill="from"), stat="identity") +
    geom_bar(data=df[I,], aes_string(x="from", y="-Freq", fill="to"), stat="identity") +
    facet_wrap(~layer) +
    scale_fill_manual(name = "Legend", values = class.colors) + 
    scale_y_continuous(labels = percent) + 
    xlab("") + 
    geom_hline(yintercept = 0) +
    coord_flip() + 
    theme(legend.position = "bottom") + 
    ylab("Percentage of land changes")
  
  gp
  
  
}



