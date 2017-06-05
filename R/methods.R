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

setGeneric("projecttwdtwRaster", 
           function(x, ...) standardGeneric("projecttwdtwRaster"))

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

writeRaster.twdtwRaster = function(x, filepath, ...){
  lapply(names(x@timeseries), function(i) writeRaster(x@timeseries[[i]], filename = paste0(filepath, "/", i, ".grd"), ...))
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
setMethod("writeRaster", "twdtwRaster", 
          definition = function(x, filepath = ".", ...) {
            writeRaster.twdtwRaster(x, filepath, ...)
          }
)

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
  if(any("doy"==layers(x)))
    return(new("twdtwRaster", timeseries=x@timeseries[i], timeline = x@timeline, doy = x@timeseries[[1]]))
  new("twdtwRaster", timeseries=x@timeseries[i], timeline = x@timeline)
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
  # if(missing(j)) j = 2:length(x@patterns)
  if(any(is.na(i))) stop("NA index not permitted")
  if(class(i)=="character") i = match(i, names(x@timeseries@timeseries))
  res = x@alignments[i]
  if(missing(j)) j = 1:length(res[[1]])
  if(class(j)=="character") j = match(j, names(x@patterns@timeseries))
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
          definition = function(object) as.character(object@labels))

#' @inheritParams twdtwTimeSeries-class
#' @rdname twdtwTimeSeries-class
#' @export
setMethod("levels", "twdtwTimeSeries",
          definition = function(x) levels(factor(labels(x))))

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
            x@timeseries = lapply(x@timeseries, crop, y=y, ...)
            x
          }
)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @param obj object of cless twdtwRaster 
#' @export
setMethod("coordinates", 
          signature = signature("twdtwRaster"),
          definition = function(obj, ...){
            coordinates(obj@timeseries[[1]], ...)
          }
)

#' @inheritParams twdtwRaster-class
#' @rdname twdtwRaster-class
#' @export
setMethod("extent", 
          signature = signature("twdtwRaster"),
          definition = function(x, y, ...){
            extent(x@timeseries[[1]])
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

# Show objects of class twdtwAssessment
show.twdtwAssessment = function(object){
  cat("An object of class \"twdtwAssessment\"\n")
  cat("Number of classification intervals:",length(object@accuracyByPeriod),"\n")
  cat("Accuracy metrics summary\n")
  cat("\nOverall\n")
  aux = object@accuracySummary$OverallAccuracy
  names(aux) = gsub("ci", "ci*", names(aux))
  print(aux, digits=2)
  cat("\nUser's\n")
  aux = object@accuracySummary$UsersAccuracy
  colnames(aux) = gsub("ci", "ci*", colnames(aux))
  print(aux, digits=2)
  cat("\nProducer's\n")
  aux = object@accuracySummary$ProducersAccuracy
  colnames(aux) = gsub("ci", "ci*", colnames(aux))
  print(aux, digits=2)
  cat("\nArea and uncertainty\n")
  aux = object@accuracySummary$AreaUncertainty
  colnames(aux) = gsub("ci", "ci*", colnames(aux))
  print(aux, digits=2)
  cat("\n*",100*object@accuracySummary$conf.int,"% confidence interval\n")
  invisible(NULL)
}

# Show objects of class twdtwCrossValidation
show.twdtwCrossValidation = function(object){
  res = summary(object, conf.int=.95)
  cat("An object of class \"twdtwCrossValidation\"\n")
  cat("Number of data partitions:",length(object@partitions),"\n")
  cat("Accuracy metrics using bootstrap simulation (CI .95)\n")
  cat("\nOverall\n")
  print(res$Overall, digits=2)
  cat("\nUser's\n")
  print(res$Users, digits=2)
  cat("\nProducer's\n")
  print(res$Producers, digits=2)
  invisible(NULL)
}

# Project raster which belong to a twdtwRaster object 
projecttwdtwRaster.twdtwRaster = function(x, to, ...){
  x@timeseries = lapply(x@timeseries, projectRaster, to, ...)
  x
}

summary.twdtwCrossValidation = function(object, conf.int=.95, ...){
  
  ov = do.call("rbind", lapply(object@accuracy, function(x){
    data.frame(OV=x$OverallAccuracy, row.names = NULL)
  }))
  
  uapa = do.call("rbind", lapply(object@accuracy, function(x){
    data.frame(label=names(x$UsersAccuracy), UA=x$UsersAccuracy, PA=x$ProducersAccuracy, row.names = NULL)
  }))
  
  sd_ov = sd(ov[, c("OV")])
  sd_uapa = aggregate(uapa[, c("UA","PA")], list(uapa$label), sd)
  l_names = levels(uapa$label)
  names(l_names) = l_names
  ic_ov = mean_cl_boot(x = ov[, c("OV")], conf.int = conf.int, ...)
  names(ic_ov) = NULL
  assess_ov = unlist(c(Accuracy=ic_ov[1], sd=sd_ov, CImin=ic_ov[2], CImax=ic_ov[3]))
  ic_ua = t(sapply(l_names, function(i) mean_cl_boot(x = uapa$UA[uapa$label==i], conf.int = conf.int, ...)))
  names(ic_ua) = NULL
  assess_ua = data.frame(Accuracy=unlist(ic_ua[,1]), sd=sd_uapa[,"UA"], CImin=unlist(ic_ua[,2]), CImax=unlist(ic_ua[,3]))
  ic_pa = t(sapply(l_names, function(i) mean_cl_boot(x = uapa$PA[uapa$label==i], conf.int = conf.int, ...)))
  names(ic_pa) = NULL  
  assess_pa = data.frame(Accuracy=unlist(ic_pa[,1]), sd=sd_uapa[,"PA"], CImin=unlist(ic_pa[,2]), CImax=unlist(ic_pa[,3]))
  list(Overall=assess_ov, Users=assess_ua, Producers=assess_pa)
}

#' @inheritParams twdtwCrossValidation-class
#' @rdname twdtwCrossValidation-class
#' @export
setMethod(f = "show", "twdtwCrossValidation",
          definition = show.twdtwCrossValidation)

#' @inheritParams twdtwAssessment-class
#' @rdname twdtwAssessment-class
#' @export
setMethod(f = "show", "twdtwAssessment",
          definition = show.twdtwAssessment)

#' @inheritParams twdtwCrossValidation-class
#' @rdname twdtwCrossValidation-class
#' @export
setMethod(f = "summary", "twdtwCrossValidation",
          definition = summary.twdtwCrossValidation)

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

#' @aliases projecttwdtwRaster
#' @inheritParams twdtwRaster-class
#' @describeIn twdtwRaster project twdtwRaster object.
#' @param crs character or object of class 'CRS'. PROJ.4 description of 
#' the coordinate reference system. For other arguments and more details see 
#' \code{\link[raster]{projectRaster}}.
#' 
#' @export
setMethod("projecttwdtwRaster", "twdtwRaster", 
          function(x, crs, ...) projecttwdtwRaster.twdtwRaster(x, crs, ...))






