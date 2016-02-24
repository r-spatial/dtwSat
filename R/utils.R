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
#   R Package dtwSat - 2016-02-22                             #
#                                                             #
###############################################################

as.list.twdtwTimeSeries = function(x) lapply(seq_along(x), function(i) x[i] )

as.list.twdtwRaster = function(x) {
    I = names(x)
    names(I) = I
    lapply(I, function(i) x[[i]])
  }

as.list.twdtwMatches = function(x) lapply(seq_along(x@alignments), 
            function(i) new("twdtwMatches", timeseries=x@timeseries[i], patterns=x@patterns, alignments=x[i, drop=FALSE]) )

dim.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = data.frame(as.character(labels(x)), t(sapply(timeseries, dim)))
  names(res) = c("label", "nrow", "ncol")
  res
}

dim.twdtwRaster = function(x){
  res = c(nlayers(x), dim=dim(x@timeseries[[1]]))
  names(res) = c("nlayers", "nrow", "ncol", "ntime")
  res
}

res.twdtwRaster = function(x){
  res(x@timeseries[[1]])
}

extent.twdtwRaster = function(x){
  extent(x@timeseries[[1]])
}

projection.twdtwRaster = function(x){
  projection(x@timeseries[[1]])
}

ncol.twdtwRaster = function(x){
  ncol(x@timeseries[[1]])
}

nrow.twdtwRaster = function(x){
  nrow(x@timeseries[[1]])
}

nlayers.twdtwRaster = function(x){
  length(names(x))
}

levels.twdtwRaster = function(x){
  levels(object@labels)
}

names.twdtwRaster = function(x){
  x@layers
}

length.twdtwRaster = function(x){
  nlayers(x)
}

index.twdtwRaster = function(x){
  x@timeline 
}

length.twdtwTimeSeries = function(x){
  length(labels(x))
}

length.twdtwMatches = function(x){
  if(length(x@alignments)<1) return(x@alignments)
  sum(sapply(x[], nrow))
}

nrow.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = sapply(timeseries, nrow)
  names(res) = as.character(labels(x))
  res
}

ncol.twdtwTimeSeries = function(x){
  timeseries = getTimeSeries(x)
  res = sapply(timeseries, ncol)
  names(res) = as.character(labels(x))
  res
}


#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "dim", "twdtwTimeSeries",
          definition = dim.twdtwTimeSeries)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "dim", "twdtwRaster",
          definition = dim.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "res", "twdtwRaster",
          definition = res.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "extent", "twdtwRaster",
          definition = extent.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "projection", "twdtwRaster",
          definition = projection.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "ncol", "twdtwRaster",
          definition = ncol.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "nrow", "twdtwRaster",
          definition = nrow.twdtwRaster)
          
#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "nlayers", "twdtwRaster",
          definition = nlayers.twdtwRaster)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "names", "twdtwRaster",
          definition = names.twdtwRaster)
          
#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "index", "twdtwRaster",
          definition = index.twdtwRaster)
          
#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "nrow", "twdtwTimeSeries",
          definition = nrow.twdtwTimeSeries)
          
#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "ncol", "twdtwTimeSeries",
          definition = ncol.twdtwTimeSeries)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "length", signature = signature("twdtwRaster"), 
          definition = length.twdtwRaster)
          
#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "length", signature = signature("twdtwTimeSeries"), 
          definition = length.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod(f = "length", signature = signature("twdtwMatches"), 
          definition = length.twdtwMatches)

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("as.list", "twdtwTimeSeries", as.list.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod("as.list", "twdtwMatches", as.list.twdtwMatches)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod("as.list", "twdtwRaster", as.list.twdtwRaster)

#' @inheritParams twdtwTimeSeries-class
#' @param i indices of the time series.
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("[", "twdtwTimeSeries", function(x, i) {
  if(missing(i)) i = 1:length(x)
  if(any(is.na(i))) stop("NA index not permitted")
  new("twdtwTimeSeries", timeseries=x@timeseries[i], labels=as.character(labels(x))[i])
})

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("[[", "twdtwTimeSeries", function(x, i) {
  if(any(is.na(i))) stop("NA index not permitted")
  x@timeseries[[i]]
})

#' @inheritParams twdtwRaster-class
#' @param i indices of the time series.
#' @rdname twdtwRaster-class
#' @export
setMethod("[", "twdtwRaster", function(x, i) {
  if(missing(i)) i = 2:nlayers(x)
  i = i[i>1]
  if(any(i > nlayers(x)))
    stop("subscript out of bounds")
  if(any(is.na(i))) stop("NA index not permitted")
  new("twdtwRaster", timeseries=x@timeseries[i], timeline = x@timeline, doy = x@timeseries[[1]])
})

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod("[[", "twdtwRaster", function(x, i) {
  if(any(is.na(i))) stop("NA index not permitted")
  x@timeseries[[i]]
})

#' @inheritParams twdtwMatches-class
#' @param i indices of the alignments.
#' @param drop if TRUE returns a data.frame, if FALSE returns a list. 
#' Default is TRUE.
#' @rdname twdtwMatches-class 
#' @export
setMethod("[", "twdtwMatches", function(x, i, drop=TRUE) {
  if(length(x@alignments)<1) return(x@alignments)
  if(missing(i)) i = 1:length(x@alignments)
  if(any(is.na(i))) stop("NA index not permitted")
  if(!drop) return(x@alignments[i])
  lapply(i, function(ii){
    xx = x@alignments[[ii]]
    res = do.call("rbind", lapply(seq_along(xx), function(j)
      data.frame(from=xx[[j]]$from, to=xx[[j]]$to, distance=xx[[j]]$distance, label=names(xx[j]))
    ))
  })
})

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class 
#' @export
setMethod("[[", c("twdtwMatches", "numeric"), function(x, i) {
  if(any(is.na(i))) stop("NA index not permitted")
  x@alignments[[i]]
})

