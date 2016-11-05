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
#' # Run TWDTW analysis for raster time series 
#' patt = MOD13Q1.MT.yearly.patterns
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' 
#' time_interval = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), 
#'                     by="12 month")
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' 
#' r_twdtw = twdtwApply(x=rts, y=patt, weight.fun=log_fun, breaks=time_interval, 
#'           filepath="~/test_twdtw", overwrite=TRUE, format="GTiff")
#' 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' plotChanges(r_lucc)
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
    geom_bar(data=df[I,], position = "stack", aes_string(x="to", y="Freq", fill="from"), stat="identity") +
    geom_bar(data=df[I,], position = "stack", aes_string(x="from", y="-Freq", fill="to"), stat="identity") +
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



