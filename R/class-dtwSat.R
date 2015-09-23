#' dtwSat class
#'
#' Class for Multidimensional Time-Weighted DTW
#'
#'
#'@section Slots :
#' \describe{
#' \item{\code{call}:}{Object of class \code{call} see 
#' \code{\link[base]{match.call}}}
#' \item{\code{alignments}:}{Object of class \code{list}. Each node of 
#' the list has one attribute, whose lengths are identical to the number of 
#' alignments. \code{from} the alignment starts, \code{to} the alignment ends, 
#' \code{distance} the DTW distance, \code{normalizedDistance} the 
#' normalized DTW distance}
#' \item{\code{mapping}:}{Object of class \code{list}. Each node of the list has 
#' the matching points of the query to the template time series.}
#' \item{\code{internals}:}{Object of class \code{list} similar to 
#' \code{\link[dtw]{dtw}} class from package \pkg{dtw}}
#' }
#' @name dtwSat
#' @aliases dtwSat-class
#' @exportClass dtwSat
#' @author Victor Maus
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
#' 
#' @description This function retrieves the results of the alignments in the object 
#' \link[dtwSat]{dtwSat}.
#' 
#' @param object A \link[dtwSat]{dtwSat} object 
#' @docType methods
#' @return object of class data.frame 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4)
#' getAlignments(alig)
#' @export
setGeneric("getAlignments", 
           function(object) {
             data.frame(object@alignments)
           }
)

#' @title Get matching points from dtwSat object
#' 
#' @description This function retrieves the matching points from the alignments 
#' in the object \link[dtwSat]{dtwSat}.
#' 
#' @param object A \link[dtwSat]{dtwSat} object 
#' @docType methods
#' @return object of class data.frame 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4)
#' getMatches(alig)
#' @export
setGeneric("getMatches", 
           function(object) {
             object@mapping
           }
)

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
#' alig = dtwSat(query.list[["Soybean"]], template, weight = "logistic", alpha = 0.1, beta = 50, alignments=4)
#' alig
#' @export
setGeneric("dtwSat", 
           function(query, template, ...) {
                twdtw(query, template, ...)
           }
)




