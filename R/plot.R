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


#' @title Plotting twdtw* objects 
#' @name plot 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Methods for plotting objects of class twdtw*.
#' 
#' @param x An object of class twdtw*.
#' @param type A character for the plot type: ''paths'', ''matches'', 
#' ''alignments'', ''classification'', ''cost'', ''patterns'', ''timeseries'',
#' ''maps'', ''area'', ''changes'', and ''distance''.
#' 
#' @param ... additional arguments to pass to plotting functions.
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotMatches}}, 
#' \code{\link[dtwSat]{plotClassification}},
#' \code{\link[dtwSat]{plotPatterns}}, 
#' \code{\link[dtwSat]{plotTimeSeries}},
#' \code{\link[dtwSat]{plotMaps}},
#' \code{\link[dtwSat]{plotArea}}, or
#' \code{\link[dtwSat]{plotChanges}}.
#'  
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @details
#' \describe{
#' 	\item{Plot types}{:
#'       \cr\code{paths}: Method for plotting the minimum paths in the cost matrix of TWDTW.
#'       \cr\code{matches}: Method for plotting the matching points from TWDTW analysis.
#'       \cr\code{alignments}: Method for plotting the alignments and respective TWDTW dissimilarity measures.
#'       \cr\code{classification}: Method for plotting the classification of each subinterval of the time series based on TWDTW analysis. 
#'       \cr\code{cost}: Method for plotting the internal matrices used during the TWDTW computation.
#'       \cr\code{patterns}: Method for plotting the temporal patterns.
#'       \cr\code{timeseries}: Method for plotting the temporal patterns.
#'       }
#' }
#' 
#' @export
NULL

#' @aliases plot-twdtwAssessment
#' @inheritParams plot
#' @rdname plot
#' @export
setMethod("plot", 
          signature(x = "twdtwAssessment"),
          function(x, type="area", ...){
            pt = pmatch(type, c("area","accuracy","map"))
            switch(pt,
                   plotAdjustedArea(x, ...),
                   plotAccuracy(x, ...),
                   plotMapSamples(x, ...)
            )
          }
)

#' @aliases plot-twdtwTimeSeries
#' @inheritParams plot
#' @rdname plot
#' @export
setMethod("plot", 
          signature(x = "twdtwCrossValidation"),
          function(x, type="crossvalidation", ...){
            pt = pmatch(type, c("crossvalidation"))
            switch(pt,
                   plotAccuracy(x, ...)
            )
          }
)

#' @aliases plot-twdtwTimeSeries
#' @inheritParams plot
#' @rdname plot
#' @export
setMethod("plot", 
          signature(x = "twdtwTimeSeries"),
          function(x, type="timeseries", ...){
            pt = pmatch(type,c("patterns","timeseries"))
            switch(pt,
                   plotPatterns(x, ...),
                   plotTimeSeries(x, ...)
            )
          }
)


#' @aliases plot-twdtwMatches
#' @inheritParams plot
#' @rdname plot
#' @export
setMethod("plot", 
          signature(x = "twdtwMatches"),
          function(x, type="alignments", ...){
            pt = pmatch(type,c("paths","matches","alignments","classification","cost"))
            switch(pt,
                   plotPaths(x, ...),
                   plotMatches(x, ...),
                   plotAlignments(x, ...),
                   plotClassification(x, ...),
                   plotCostMatrix(x, ...)
            )
          }
)



#' @aliases plot-twdtwRaster
#' @inheritParams plot
#' @rdname plot
#' @export
setMethod("plot", signature(x = "twdtwRaster"), function(x, type="maps", ...) .PlotRaster(x, type=type, ...))

.PlotRaster = function(x, type, time.levels=NULL, time.labels=NULL, class.levels=NULL, class.labels=NULL, class.colors=NULL, layers=NULL, perc=TRUE, ...){
  
  if(type=="distance") {
      
      if( is.null(time.levels))
         time.levels = names(x)
          
      if(is(time.levels, "numeric"))
         time.levels = names(x)[time.levels]

      if( is.null(time.labels))
         time.labels = format(as.Date(time.levels, "date.%Y.%m.%d"), "%Y")
            
      if(length(time.levels)!=length(time.labels))
        stop("time.levels and time.labels have different length")
        
      if(is.null(layers)) {
        if(any(coverages(x)=="Distance")){
          layers = "Distance" 
          time.levels = time.levels 
          labels = time.labels 
          time.labels=NULL
        }
        else {
          layers = coverages(x)
          layers = layers[!layers%in%"doy"]
          time.levels = time.levels[1]
          time.labels = time.labels[1]
          labels = layers
        }
      } else {
          if(is(layers, "numeric")) layers = coverages(x)[layers]
          time.levels = time.levels[1]
          time.labels = time.labels[1]
          labels = layers
      }
      
      x = lapply(as.list(x)[layers], FUN=subset, subset=time.levels)
      gp = .plotDistance(brick(x), layers, labels, time.labels)
      
  } else {
    
      if( is.null(time.levels))
        time.levels = seq_along(index(x))
      
      if(is.null(time.labels))
        time.labels = format(index(x), "%Y")
      
      if(is(time.levels, "numeric")){
        time.levels = names(x)[time.levels]
        time.labels = time.labels[time.levels]
      }
        
      
      if(length(time.levels)!=length(time.labels))
        stop("time.levels and time.labels have different length")
      
      # if(length(time.levels)>16){
      #   time.levels = time.levels[1:16]
      #   time.labels = time.labels[1:16]
      # }
      
      if( is.null(class.levels))
        class.levels = levels(x)
        
      if(length(class.levels)<1)
        class.levels = sort(unique(as.numeric(x[["Class"]][])))
      
      if( is.null(class.labels))  
        class.labels = labels(x)
          
      if(length(class.labels)<1)
        class.labels = as.character(class.levels)
      
      if( is.null(class.colors) )
        class.colors = brewer.pal(length(class.levels), "Set3")
      
      if( length(class.colors)<length(class.levels) )
        class.colors = brewer.pal(length(class.levels), "Set3")
      
      if(length(class.levels)!=length(class.labels))
        stop("class.levels and class.labels have different length")  

      names(time.labels)  = time.levels
      names(time.levels)  = time.labels
      names(class.colors) = class.labels
      names(class.levels) = class.labels
      names(class.labels) = class.labels
      
      x = subset(x=x[["Class"]], subset=time.levels)

      pt = pmatch(type, c("maps","area","changes"))
      
      gp = switch(pt,
            .plotMaps(x, time.levels, time.labels, class.levels, class.labels, class.colors),
            .plotArea(x, time.levels, time.labels, class.levels, class.labels, class.colors, perc),
            .plotChanges(x, time.levels, time.labels, class.levels, class.labels, class.colors)
      )
  }
  gp
}


