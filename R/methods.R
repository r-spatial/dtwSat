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


setGeneric("layers", 
           function(x) standardGeneric("layers"))

setGeneric("coverages", 
           function(x) standardGeneric("coverages"))

setGeneric("bands", 
           function(x) standardGeneric("bands"))
           
setGeneric("is.twdtwTimeSeries", 
           function(x) standardGeneric("is.twdtwTimeSeries"))

setGeneric("is.twdtwMatches", 
           function(x) standardGeneric("is.twdtwMatches"))

setGeneric("is.twdtwRaster", 
           function(x) standardGeneric("is.twdtwRaster"))
           
as.list.twdtwTimeSeries = function(x) lapply(seq_along(x), function(i) 
      new("twdtwTimeSeries", x[[i]], labels(x)[i]) )

as.list.twdtwRaster = function(x) {
    I = coverages(x)
    names(I) = I
    lapply(I, function(i) x[[i]])
  }

as.list.twdtwMatches = function(x) lapply(seq_along(x@timeseries), function(i) 
      new("twdtwMatches", new("twdtwTimeSeries", x@timeseries[[i]], labels(x@timeseries)[i]), x@patterns, list(x@alignments[[i]])) )

dim.twdtwTimeSeries = function(x){
  res = data.frame(as.character(labels(x)), t(sapply(x@timeseries, dim)))
  names(res) = c("label", "nrow", "ncol")
  row.names(res) = NULL
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
  length(coverages(x))
}

levels.twdtwRaster = function(x){
  x@levels
}

layers.twdtwRaster = function(x){
  x@layers
}

coverages.twdtwRaster = function(x){
  x@layers
}

bands.twdtwRaster = function(x){
  x@layers
}

names.twdtwRaster = function(x){
  names(x@timeline)
}

length.twdtwRaster = function(x){
  nlayers(x)
}

index.twdtwRaster = function(x){
  x@timeline 
}

index.twdtwTimeSeries = function(x){
  lapply(x@timeseries, index)
}

index.twdtwMatches = function(x){
  lapply(getTimeSeries(x), index)
}

length.twdtwTimeSeries = function(x){
  length(x@timeseries)
}

length.twdtwMatches = function(x){
  if(length(x@alignments)<1) return(x@alignments)
  sum(sapply(x[], nrow))
}

nrow.twdtwTimeSeries = function(x){
  res = sapply(x@timeseries, nrow)
  names(res) = as.character(labels(x))
  res
}

ncol.twdtwTimeSeries = function(x){
  res = sapply(x@timeseries, ncol)
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
setMethod(f = "levels", "twdtwRaster",
          definition = levels.twdtwRaster)

#' @aliases layers
#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "layers", "twdtwRaster",
          definition = layers.twdtwRaster)

#' @aliases coverages
#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "coverages", "twdtwRaster",
          definition = coverages.twdtwRaster)
          
#' @aliases bands
#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "bands", "twdtwRaster",
          definition = bands.twdtwRaster)
         
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
setMethod(f = "index", "twdtwTimeSeries",
          definition = index.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod(f = "index", "twdtwMatches",
          definition = index.twdtwMatches)
          
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
  x@timeseries[i]
})

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("[[", "twdtwTimeSeries", function(x, i) {
  if(any(is.na(i))) stop("NA index not permitted")
  x@timeseries[[i, drop=FALSE]]
})

#' @inheritParams twdtwRaster-class
#' @param i indices of the time series.
#' @rdname twdtwRaster-class
#' @export
setMethod("[", "twdtwRaster", function(x, i) {
  if(missing(i)) i = 2:nlayers(x)
  if(length(i)>1) i = i[i>1]
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
#' @param i indices of the time series.
#' @param j indices of the pattern.
#' @param drop if TRUE returns a data.frame, if FALSE returns a list. 
#' Default is TRUE.
#' @rdname twdtwMatches-class 
#' @export
setMethod("[", "twdtwMatches", function(x, i, j, drop=TRUE) {
  if(length(x@alignments)<1) return(x@alignments)
  if(missing(i)) i = 1:length(x@alignments)
  if(any(is.na(i))) stop("NA index not permitted")
  res = x@alignments[i]
  if(missing(j)) j = 1:length(res[[1]])
  if(any(is.na(j))) stop("NA index not permitted")
  res = lapply(res, function(x) x[j])
  res = res[sapply(res, length)>0]
  if(!drop) return(res)
  lapply(res, function(x){
    res = do.call("rbind", lapply(seq_along(x), function(jj){
            data.frame(Alig.N=seq_along(x[[jj]]$distance),from=x[[jj]]$from, to=x[[jj]]$to, distance=x[[jj]]$distance, label=x[[jj]]$label, row.names=NULL)
          }))
    res[order(res$from),]      
  })
})

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class 
#' @export
setMethod("[[", c("twdtwMatches", "numeric"), function(x, i, j,drop=TRUE) {
  if(any(is.na(i))) stop("NA index not permitted")
  if(missing(j)) j = 1:length(x@alignments[[1]])
  x[i,j,drop=drop][[1]]
})

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("labels", signature = signature(object="twdtwTimeSeries"),
          definition = function(object) object@labels)

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("levels", "twdtwTimeSeries",
          definition = function(x) levels(labels(x)))

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod("labels", signature = signature(object="twdtwRaster"),
          definition = function(object) as.character(object@labels))
          
#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod("labels", 
          signature = signature(object="twdtwMatches"),
          definition = function(object){
            list(timeseries = labels(object@timeseries), 
                 patterns = labels(object@patterns))
          }
)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @param y Extent object, or any object from which an Extent object can be extracted.
#' @export
setMethod("crop", 
          signature = signature("twdtwRaster"),
          definition = function(x, y, ...){
            x@timeseries = lapply(x@timeseries, crop, y=extent(y), ...)
            x
          }
)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod("extent", 
          signature = signature("twdtwRaster"),
          definition = function(x, y, ...){
            extent(x@timeseries$doy)
          }
)

# Show objects of class twdtwTimeSeries 
show.twdtwTimeSeries = function(object){
  cat("An object of class \"twdtwTimeSeries\"\n")
  cat("Slot \"timeseries\" length:",length(object),"\n")
  cat("Slot \"labels\": ")
  I = match(1:3, seq_along(labels(object)))
  print(labels(object)[na.omit(I)])
  invisible(NULL)
}

# Show objects of class twdtwMatches 
show.twdtwMatches = function(object){
  cat("An object of class \"twdtwMatches\"\n")
  cat("Number of time series:",length(object@timeseries),"\n")
  cat("Number of Alignments:",length(object),"\n")
  cat("Patterns labels:",as.character(labels(object@patterns)),"\n")
  invisible(NULL)
}

# Show objects of class twdtwRaster 
show.twdtwRaster = function(object){
  cat("An object of class \"twdtwRaster\"\n")
  cat("Time series layers:",coverages(object),"\n")
  cat("Time range:",paste(min(object@timeline)),"...",paste(max(object@timeline)),"\n")
  cat("dimensions:",dim(object),"(nlayers, nrow, ncol, length)\n")
  cat("resolution:",res(object)," (x, y)\n")
  cat("extent    :",as.vector(extent(object)), "(xmin, xmax, ymin, ymax)\n")
  cat("coord.ref.:",projection(object),"\n") 
  invisible(NULL)
}

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod(f = "show", "twdtwTimeSeries",
          definition = show.twdtwTimeSeries)

#' @inheritParams twdtwMatches-class
#' @rdname twdtwMatches-class
#' @export
setMethod(f = "show", "twdtwMatches",
          definition = show.twdtwMatches)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod(f = "show", "twdtwRaster",
          definition = show.twdtwRaster)
           
#' @aliases is.twdtwTimeSeries
#' @inheritParams twdtwTimeSeries-class
#' @describeIn twdtwTimeSeries Check if the object belongs to the class twdtwTimeSeries.
#' @export
setMethod("is.twdtwTimeSeries", "ANY", 
          function(x) is(x, "twdtwTimeSeries"))

#' @aliases is.twdtwMatches
#' @inheritParams twdtwMatches-class
#' @describeIn twdtwMatches Check if the object belongs to the class twdtwMatches.
#' @export
setMethod("is.twdtwMatches", "ANY", 
          function(x) is(x, "twdtwMatches"))
          
#' @aliases is.twdtwRaster
#' @inheritParams twdtwRaster-class
#' @describeIn twdtwRaster Check if the object belongs to the class twdtwRaster.
#' @export
setMethod("is.twdtwRaster", "ANY", 
          function(x) is(x, "twdtwRaster"))
         
# #' @aliases summary
# #' @inheritParams twdtwMatches-class
# #' @describeIn twdtwMatches Summary of objects of class twdtwMatches.
# #' @export         
# setMethod("summary", 
#           signature(object = "twdtwMatches"),
#           function(object, labels=NULL, ...){
#             summary.twdtw(object, labels=labels, ...)
#           }
# )
# 
# summary.twdtw = function(object, ...){
#   lapply(as.list(object), function(obj){
#     res = lapply(labels(obj)$patterns, function(l){
#       m = subset(obj, patterns.labels=l)[[1]]
#       c(labels=l, N.Matches=nrow(m), summary(m$distance))
#     })
#     data.frame(res)
#   })
# }

