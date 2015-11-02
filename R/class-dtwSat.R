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
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### dtwSat CLASS AND METHODS
#' @title dtwSat-class
#' @name dtwSat-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for Multidimensional Time-Weighted DTW alignments
#' 
#' 
#' @section Slots :
#' \describe{
#'  \item{\code{call}:}{An object of class \code{\link[base]{call}}, see 
#'  \code{\link[base]{match.call}}.}
#'  \item{\code{alignments}:}{A named \code{\link[base]{list}} whose elements 
#'  have length identical to the number of alignments.
#'  The elements are:  
#'       \cr\code{query}: query identification,
#'       \cr\code{from}: starting dates,
#'       \cr\code{to}: ending dates, and
#' 	     \cr\code{distance}: TWDTW distances.
#'       }
#'  \item{\code{mapping}:}{An object of class \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the query and the time series. 
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the query, and
#'       \cr\code{index2}: matching points of the time series.
#'       }
#' 	\item{\code{internals}:}{An object of class \code{\link[base]{list}} whose 
#'   elements have the internal structures used by \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken in the alignments,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}}
#'       \cr\code{query}: query time series, 
#'       \cr\code{timeseries}: satellite image time series,
#'       \cr\code{N}: \code{query} length, and 
#'       \cr\code{M}: \code{timeseries} length.
#'       }
#' }
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtw]{stepPattern}}, and \code{\link[dtw]{dtw}}.
#' 
NULL
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
                                query=numeric(0), 
                                from=numeric(0), 
                                to=numeric(0), 
                                distance=numeric(0))
      .Object@mapping =    list(
                                index1 = numeric(0), 
                                index2 = numeric(0))
      .Object@internals =  list(
                                timeWeight = matrix(0),
                                localMatrix = matrix(0),
                                costMatrix = matrix(0),
                                directionMatrix = matrix(0),
                                stepPattern = numeric(0),
                                N = numeric(0),
                                M = numeric(0),
                                query = numeric(0),
                                timeseries = numeric(0)
                                )
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

setMethod("show", 
          signature = signature(object="dtwSat"),
          definition = function(object){
            cat("Time-Weighted DTW alignment object\n")
            cat("Number of alignments:",nrow(getAlignments(object)),"\n")
            print(head(getAlignments(object)))
            invisible(NULL)
          }
)




###############################################################
#### GENERIC METHODS


#' @title Multidimensional Time-Weighted DTW Alignment
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves one or more possible alignments of a query within 
#' a time series.
#' 
#' @param query A \link[zoo]{zoo} object with a query time series.
#' @param timeseries A \link[zoo]{zoo} object with a time series similar 
#' to \code{query}. \code{timeseries} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' \emph{i.e.} \code{nrow(query)<=nrow(timeseries)}.
#' @param ... additional arguments passed to \code{\link[dtwSat]{twdtw}}
#' @param template is deprecated, please use \code{timeseries} instead
#' 
#' @docType methods
#' @return An object of class \link[dtwSat]{dtwSat-class}
#'  
#' @seealso \code{\link[dtwSat]{twdtw}} and \code{\link[dtwSat]{mtwdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template)
#' alig
#' @export
#' @name dtwSat
setGeneric("dtwSat", 
           function(query, timeseries, template, ...) {
             if (!missing(template)){
               warning("argument template is deprecated, please use timeseries instead", call. = FALSE)
               timeseries = template
             }
             twdtw(query, timeseries, ...)
           }
)


#' @title Get alignments from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the alignments 
#' in the object \link[dtwSat]{dtwSat-class}
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @docType methods
#' @return An object of class \code{\link[base]{list}} whose 
#' elements have length identical to the number of alignments.
#' The elements are:  
#'       \cr\code{query}: query identification,
#'       \cr\code{from}: starting dates,
#'       \cr\code{to}: ending dates, and
#' 	     \cr\code{distance}: TWDTW distances.
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}} and \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template)
#' getAlignments(alig)
#' @export
setGeneric("getAlignments", 
           function(object) {
             data.frame(object@alignments, stringsAsFactors = FALSE)
           }
)

#' @title Get matching points from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the matching points 
#' for each alignment between the \code{query} and the \code{timeseries}
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @docType methods
#' @return An object of class \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the query and the time series. 
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the query, and
#'       \cr\code{index2}: matching points of the time series.
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template)
#' getMatches(alig)
#' @export
setGeneric("getMatches", 
           function(object) {
             object@mapping
           }
)



#' @title Get internals from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves cost matrix, inputs, and other 
#' internal structures from dtwSat-class object
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @docType methods
#' @return An object of class \code{\link[base]{list}} whose 
#'   elements have the internal structures used by \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken in the alignments,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}}
#'       \cr\code{query}: query time series, 
#'       \cr\code{timeseries}: satellite image time series,
#'       \cr\code{N}: \code{query} length, and 
#'       \cr\code{M}: \code{timeseries} length.
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], timeseries=template, keep=TRUE)
#' getInternals(alig)
#' @export
setGeneric("getInternals", 
           function(object){
             object@internals
           }
)


