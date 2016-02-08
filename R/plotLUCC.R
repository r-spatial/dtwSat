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


#' @title Plotting land use/cover changes (LUCC) 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting the land use and land cover 
#' changes 
#' 
#' @param x A \code{\link[raster]{Raster-class}}
#' \code{\link[raster]{brick}} or \code{\link[raster]{stack}} object. 
#' A \link[base]{character} with the raster file name is also accepted.
#' @param type A \link[base]{character}. The plot type: ''map'' for raster maps, 
#' ''area'' for the accumulated area over time, or ''change'' for the land changes.
#' Default is ''map''.
#' @param layer.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the layers to plot. For plot type ''change'' the minimum length 
#' is two.
#' @param layer.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the layers. It must have the same 
#' length as layer.levels. Default is NULL.
#' @param class.levels A \link[base]{character} or \link[base]{numeric}
#' vector with the levels of the raster values. Default is NULL. 
#' @param class.labels A \link[base]{character} or \link[base]{numeric}
#' vector with the labels of the raster values. It must have the same 
#' length as class.levels. Default is NULL.
#' @param class.colors a set of aesthetic values. It must have the same 
#' length as class.levels. Default is NULL. See 
#' \link[ggplot2]{scale_fill_manual} for details.
#' 
#' @docType methods
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwApply}}, 
#' \code{\link[raster]{brick}}, and
#' \code{\link[raster]{stack}}.
#'  
#' @examples
#' 
#' levels = c(seq_along(patterns_vignette.list), 255)
#' labels = c(names(patterns_vignette.list), "Unclassified")
#' colors = c("#996400", "#005500", "#D8B777", "#E6D219", "#E6BEC8", "#C8C8C8")
#' names(colors) = labels
#' 
#' x = system.file('lucc_MT.tif',  package = 'dtwSat')
#' 
#' ### Plot maps
#' #gp1 = plotLUCC(x = x, type = "map", layer.labels = 2008:2013, 
#' #        class.levels = levels, class.labels = labels, 
#' #        class.colors = colors)
#' #gp1 
#' 
#' # Plot area 
#' #gp2 = plotLUCC(x = x, type = "area", layer.labels = 2008:2013, 
#' #        class.levels = levels, class.labels = labels, 
#' #        class.colors = colors)
#' #gp2
#' 
#' # Plot land use changes 
#' #gp3 = plotLUCC(x = x, type = "change", layer.labels = 2008:2013, 
#' #        class.levels = levels, class.labels = labels, 
#' #        class.colors = colors)
#' #gp3
#' 
#' @export
plotLUCC = function(x, type="area", layer.levels, layer.labels=NULL, 
                    class.levels=NULL, class.labels=NULL, class.colors=NULL){
  
  if(is(x, "character"))
    x = brick(x)
  
  if(!is(x,"RasterBrick") & !is(x,"RasterStack"))
    stop("x is not a RasterBrick or RasterStack object")
  
  if(missing(layer.levels))
    layer.levels = names(x)
  
  if(is.null(layer.labels))
    layer.labels = layer.levels
  
  if(length(layer.levels)!=length(layer.labels))
    stop("layer.levels and layer.labels have different length")
  
  x_levels = unique(as.vector(unique(x[])))
  if( is.null(class.colors) | length(class.colors)<length(x_levels) )
    class.colors = brewer.pal(length(x_levels), "Set3")
  
  if(is.null(class.levels))
    class.levels = sort(x_levels)
  
  if(is.null(class.labels))
    class.labels = class.levels
  
  if(length(class.levels)!=length(class.labels))
    stop("class.levels and class.labels have different length")  
  
  names(class.colors) = class.labels
  names(class.levels) = class.labels
  
  pt = pmatch(type,c("map","area","change"))
  
  x = raster::subset(x=x, subset=layer.levels)
  names(layer.levels) = names(x)
  names(layer.labels) = names(x)
  
  switch(pt,
         .plotLUMap(x, layer.levels, layer.labels, class.levels, class.labels, class.colors),
         .plotLUArea(x, layer.levels, layer.labels, class.levels, class.labels, class.colors),
         .plotLUChange(x, layer.levels, layer.labels, class.levels, class.labels, class.colors)
  )
  
}

.plotLUChange = function(x, layer.levels, layer.labels, class.levels, class.labels, class.colors){
  
  if(length(layer.levels)<2)
    stop("the vector layer.levels is shorter than two")
  
  df = do.call("rbind", lapply(seq_along(layer.levels)[-1], function(l){
    from = raster::subset(x=x, subset=layer.levels[l-1])[]
    to   = raster::subset(x=x, subset=layer.levels[l]  )[]
    res = data.frame(from, to)
    res = data.frame(xtabs(~ from + to, res) / nrow(res))
    res$layer = paste0(layer.labels[l-1],"-",layer.labels[l])
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

.plotLUMap = function(x, layer.levels, layer.labels, class.levels, class.labels, class.colors){
  
  df.map = data.frame(coordinates(x), x[] )
  df.map = melt(df.map, id.vars = c("x", "y"))
  df.map$value = factor(df.map$value, levels = class.levels, labels = class.labels)
  
  df.map$variable = factor(df.map$variable, labels = layer.labels)
  
  gp = ggplot(data=df.map, aes_string(x="x", y="y")) +
    geom_raster(aes_string(fill="value")) + 
    scale_fill_manual(name="Legend", values = class.colors) + 
    facet_wrap(~variable) + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) + 
    theme(legend.position = "bottom") + 
    coord_fixed(ratio = 1) + 
    xlab("Longitude") + 
    ylab("Latitude")
  gp 
  
}

.plotLUArea = function(x, layer.levels, layer.labels, class.levels, class.labels, class.colors){
  
  df.map = data.frame(coordinates(x), x[] )
  df.map = melt(df.map, id.vars = c("x", "y"))
  df.map$value = factor(df.map$value, levels = class.levels, labels = class.labels)
  
  df.area = data.frame(prop.table(xtabs(~ variable + value, df.map), 1))
  # df.area$Time = factor(as.numeric(df.area$variable), labels = layer.labels)
  df.area$Time = as.numeric(df.area$variable)
  
  x.breaks = pretty_breaks()(range(df.area$Time, na.rm = TRUE))
  
  gp = ggplot(data=df.area, aes_string(x="Time", y="Freq", fill="value")) +
    geom_area(position = 'stack') + 
    scale_fill_manual(name="Legend", values = class.colors) + 
    scale_x_continuous(expand = c(0.01, 0), breaks = x.breaks, labels = layer.labels) + 
    scale_y_continuous(expand = c(0, 0), labels = percent) +
    theme(legend.position = "bottom",
          panel.background = element_blank()) + 
    ylab("Percentage of total land area")
  gp 
  
}

