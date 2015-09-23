#' @title dtwSat-class
#'   
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for Multidimensional Time-Weighted DTW
#' 
#' @import methods
#' @import dtw
#' @import zoo
#' @importFrom proxy dist
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @importFrom waveslim mra
#' @importFrom ggplot2 ggplot geom_line geom_point geom_path geom_raster xlab ylab scale_x_continuous scale_y_continuous scale_x_date scale_y_date annotate scale_fill_gradientn aes_string waiver
#' @importFrom scales pretty_breaks
#' @useDynLib dtwSat computeCM
#'
#' @section Slots :
#' \describe{
#' 	\item{\code{call}:}{Object of class \code{call} see 
#' 	\code{\link[base]{match.call}}.}
#' 	\item{\code{alignments}:}{Object of class \code{list} where each node of 
#' 	the list has one attribute, whose lengths are identical to the number of 
#' 	alignments. The attributes are:} 
#' 	  \item{\code{from}}{ the alignment starts, }
#' 	  \item{\code{to}}{ the alignment ends, }
#' 	  \item{\code{distance}}{ the DTW distance, and}
#' 	  \item{\code{normalizedDistance}}{ the normalized DTW distance.}
#' 	\item{\code{mapping}:}{ Object of class \code{list}. Each node of the list has 
#' 	the matching points of the query to the template time series.}
#' 	\item{\code{internals}:}{ Object of class \code{list} similar to 
#' 	\code{\link[dtw]{dtw}} class from package \pkg{dtw}.}
#' }
#' @aliases dtwSat-class
dtwSat = setClass(
  Class = "dtwSat",
  slots = c(
    call = "call",
    alignments = "list",
    mapping = "list",
    internals = "list"
  ),
  validity = function(object){
    if(!is(object@call, "call")){
      stop("[dtwSat: validation] Invalid call object, class different from call.")
    }else{}
    if(!is(object@alignments, "list")){
      stop("[dtwSat: validation] Invalid TWDTW alignments, class different from data.frema.")
    }else{}
    if(!is(object@mapping, "list")){
      stop("[dtwSat: validation] Invalid TWDTW mapping, class different from list.")
    }else{}
    if(!is(object@internals, "list")){
      stop("[dtwSat: validation] Invalid DTW internals, class different from dtw.")
    }else{}
    return(TRUE)
  }
)


#' @title dtwSat-initialize
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Initialize object of class dtwSat
#'
#' @description Show object of class dtwSat
#' @param .Object A character "dtwSat"
#' @param call Object of class \code{call} see \code{\link[base]{match.call}}
#' @param internals Object of class \code{list}
#' @param alignments Object of class \code{list}
#' @param mapping Object of class \code{list}
#' 
#' @docType methods
#' 
setMethod("initialize",
  signature = "dtwSat",
  definition = 
    function(.Object, call, internals, alignments, mapping){
      .Object@call = new("call")
      .Object@alignments = list(
                                from=numeric(0), 
                                to=numeric(0), 
                                distance=numeric(0), 
                                normalizedDistance=numeric(0))
      .Object@mapping =    list(
                                index1 = numeric(0), 
                                index2 = numeric(0))
      .Object@internals =  list()
      if(!missing(call))
        .Object@call = call
      if(!missing(alignments))
        .Object@alignments = alignments
      if(!missing(mapping))
        .Object@mapping = mapping
      if(!missing(internals))
        .Object@internals = internals
      validObject(.Object)
      return(.Object)
  }
)

#' @title Get alignments from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the results of the alignments in the object 
#' \link[dtwSat]{dtwSat}
#' 
#' @param object A \link[dtwSat]{dtwSat} object
#' @docType methods
#' @return An object of class data.frame
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'              alpha = 0.1, beta = 50, alignments=4)
#' getAlignments(alig)
#' @export
setGeneric("getAlignments", 
           function(object) {
             data.frame(object@alignments)
           }
)

#' @title Get matching points from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the matching points 
#' from the alignments in the object \link[dtwSat]{dtwSat}
#' 
#' @param object A \link[dtwSat]{dtwSat} object
#' @docType methods
#' @return object of class data.frame 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 50, alignments=4)
#' getMatches(alig)
#' @export
setGeneric("getMatches", 
           function(object) {
             object@mapping
           }
)


#' @title dtwSat-show
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Show object of class dtwSat
#' @param object A \link[dtwSat]{dtwSat} object
#' @docType methods
#' 
setMethod("show", 
          signature = signature(object="dtwSat"),
          definition = function(object){
            cat("Time-Weighted DTW alignment object\n")
            cat("Alignments:\n")
            print(getAlignments(object))
            invisible(NULL)
          }
)


#' @title Multidimensional Time-Weighted DTW analysis
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves one or more possible alignments of a query within 
#' a time series.
#' 
#' @param query A \link[zoo]{zoo} object with the multidimensional time series.
#' @param template A \link[zoo]{zoo} object with the template time series. The index of 
#' the zoo object must be of class \code{\link[base]{Date}}.
#' It must be iguel or be equal or longer than the length of the query and 
#' the same number of dimensions. The index of the zoo object must be of 
#' class \code{\link[base]{Date}}.
#' @param ... see \code{\link[dtwSat]{twdtw}}
#' @docType methods
#' @return object of class \code{\link[dtwSat]{dtwSat}} 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 50, alignments=4)
#' alig
#' @export
#' @name dtwSat
setGeneric("dtwSat", 
           function(query, template, ...) {
                twdtw(query, template, ...)
           }
)




