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
#   R Package dtwSat - 2016-02-18                             #
#                                                             #
###############################################################

#' @title Subset time series 
#' @name subset
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Get subsets from objects of class twdtw*.
#' 
#' @inheritParams get 
#' @param x An objects of class twdtw*.
#' 
#' @param k A positive integer. The index of the last alignment to include in 
#' the subset.
#' 
#' @param e An extent object, or any object from which an Extent object can
#' be extracted. See \link[raster]{crop} for details.
#' 
#' @param layers a vector with the names of the \code{twdtwRaster} object to include in 
#' the subset.
#' 
#' @param labels character vector with time series labels.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwRaster-class}}, 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}, and 
#' \code{\link[dtwSat]{twdtwMatches-class}}
#'
#' @return an object of class twdtw*. 
#'
#' @examples
#' # Getting time series from objects of class twdtwTimeSeries
#' ts = twdtwTimeSeries(MOD13Q1.ts.list)
#' ts = subset(ts, 2)
#' ts
#' # Getting time series from objects of class twdtwTimeSeries
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' mat = twdtwApply(x=ts, y=patt, weight.fun=logisticWeight(-0.1,100))
#' mat = subset(mat, k=4)
#' 
#' ## This example creates a twdtwRaster object and extract time series from it. 
#'
#' # Creating objects of class twdtwRaster with evi and ndvi time series 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, timeline=timeline)
#' 
#' rts_evi = subset(rts, layers="evi")
#'
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' prj_string = scan(system.file("lucc_MT/data/samples_projection", package="dtwSat"), 
#'                   what = "character")
#' 
#' # Extract time series 
#' ts_evi = getTimeSeries(rts_evi, y = field_samples, proj4string = prj_string)
#' 
#' # subset all labels = "Forest"
#' ts_forest = subset(ts_evi, labels="Forest")
#' 
NULL

#' @aliases subset-twdtwTimeSeries
#' @inheritParams subset
#' @rdname subset
#' @export
setMethod("subset", "twdtwTimeSeries", function(x, labels=NULL) 
          subset.twdtwTimeSeries(x=x, labels=labels))


subset.twdtwTimeSeries = function(x, labels){
  if(is.null(labels)) labels = labels(x)
  if(is.numeric(labels)) return(twdtwTimeSeries(x@timeseries[labels], labels=x@labels[labels]))
  I = which(!is.na(match(x@labels, labels)))
  if(length(I)<1) return(new("twdtwTimeSeries"))
  twdtwTimeSeries(x@timeseries[I], labels=x@labels[I])
}


#' @aliases subset-twdtwMatches
#' @inheritParams subset
#' @rdname subset
#' @export
setMethod("subset", "twdtwMatches", function(x, timeseries.labels=NULL, patterns.labels=NULL, k=NULL) 
          subset.twdtwMatches(x=x, timeseries.labels=timeseries.labels, patterns.labels=patterns.labels, k=k) )


subset.twdtwMatches = function(x, timeseries.labels, patterns.labels, k){
  if(is.null(timeseries.labels)) timeseries.labels = as.character(labels(x@timeseries))
  if(is.null(patterns.labels)) patterns.labels = as.character(labels(x@patterns))
  if(is.null(k)) k = 1:length(x)
  k = unique(k)
  I = timeseries.labels
  J = patterns.labels
  if(is.character(I)) I = which(!is.na(match(x@timeseries@labels, timeseries.labels)))
  if(is.character(J)) J = which(!is.na(match(x@patterns@labels, patterns.labels)))
  timeseries = subset(x@timeseries, labels=I)
  patterns = subset(x@patterns, labels=J)
  names(J) = labels(patterns)
  alignments = lapply(I, function(i){
    out = lapply(J, function(j){
      res = x@alignments[[i]][[j]]
      k = k[ k<=res$K ]
      res$K = length(k)
      res$from = res$from[k]
      res$to = res$to[k]
      res$distance = res$distance[k]
      if(length(k)<1) res$label = numeric(0)
      if(length(res$matching)>length(k)) res$matching = res$matching[k]
      res
    })
    # names(out) = patterns.labels
    out
  })
  twdtwMatches(timeseries=timeseries, patterns=patterns, alignments=alignments)
}

#' @aliases subset-twdtwRaster
#' @inheritParams subset
#' @rdname subset
#' @export
setMethod("subset", "twdtwRaster", function(x, e=NULL, layers=NULL) 
          subset.twdtwRaster(x=x, e=e, layers=layers) )

subset.twdtwRaster = function(x, e, layers){
    if(is.null(layers)) 
      layers = names(x)
    if(is.null(e))
      e = extent(x)
    res = x
    res@layers = layers
    res@timeseries = res@timeseries[layers]
    res@timeseries = lapply(res@timeseries, crop, y=e)
    res
}






          
          
          
          
          
          
          
          
          
          
          
          
          