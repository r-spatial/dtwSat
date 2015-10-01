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
#'       \cr\code{query}: a query identification,
#'       \cr\code{from}: starting dates,
#'       \cr\code{to}: ending dates,
#' 	     \cr\code{distance}: TWDTW distances, and
#' 	     \cr\code{normalizedDistance}: normalized DTW distances.
#'        }
#'  \item{\code{mapping}:}{An object of class \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the query and the template time series. 
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the query, and
#'       \cr\code{index2}: matching points of the template.
#'       }
#' 	\item{\code{internals}:}{An object of class \code{\link[base]{list}} whose 
#'   elements have the internal structures used by \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{costMatrix}: the cumulative cost matrix,
#'       \cr\code{stepPattern}: the \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}},
#'       \cr\code{N}: query length
#'       \cr\code{M}: reference length
#'       \cr\code{query}: the query time series, and
#'       \cr\code{template}: the reference time series.
#'       }
#' }
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtw]{stepPattern}}, and \code{\link[dtw]{dtw}}.
#' 
#' @import methods
#' @import dtw
#' @import zoo
#' @importFrom proxy dist
#' @importFrom reshape2 melt
#' @importFrom graphics plot
#' @importFrom waveslim mra
#' @importFrom ggplot2 ggplot geom_line geom_point geom_path geom_raster geom_polygon xlab ylab scale_x_continuous scale_y_continuous scale_x_date scale_y_date scale_fill_brewer annotate scale_fill_gradientn aes_string waiver
#' @importFrom scales pretty_breaks
#' @useDynLib dtwSat computeCM
#' @importFrom grDevices terrain.colors
#' @importFrom utils tail head
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
                                distance=numeric(0), 
                                normalizedDistance=numeric(0))
      .Object@mapping =    list(
                                index1 = numeric(0), 
                                index2 = numeric(0))
      .Object@internals =  list(
                                costMatrix = matrix(0),
                                stepPattern = numeric(0),
                                N = numeric(0),
                                M = numeric(0),
                                query = numeric(0),
                                template = numeric(0)
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
            cat("Alignments:\n")
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
#' @param query A \link[zoo]{zoo} object with a query time series
#' @param template A \link[zoo]{zoo} object with a template time series similar 
#' to \code{query}. The \code{template} must have the same number of attributes
#' and be equal to or longer than the \code{query}
#' @param ... additional arguments passed to \code{\link[dtwSat]{twdtw}}
#' @docType methods
#' @return An object of class \link[dtwSat]{dtwSat-class}
#'  
#' @seealso \code{\link[dtwSat]{twdtw}} and \code{\link[dtwSat]{mtwdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template)
#' alig
#' @export
#' @name dtwSat
setGeneric("dtwSat", 
           function(query, template, ...) {
             twdtw(query, template, ...)
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
#' @return An object of class \link[base]{data.frame}
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}} and \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template)
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
#' for each alignment between the \code{query} and the \code{template}
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @docType methods
#' @return An object of class \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the query and the template time series.
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the query, and
#'       \cr\code{index2}: matching points of the template
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template)
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
#'  elements are internals from \link[dtwSat]{dtwSat-class} object 
#'       \cr\code{costMatrix}: the cumulative cost matrix,
#'       \cr\code{stepPattern}: the \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}},
#'       \cr\code{N}: query length
#'       \cr\code{M}: reference length
#'       \cr\code{query}: the query time series, and
#'       \cr\code{template}: the reference time series.
#' 
#' @seealso \code{\link[dtwSat]{dtwSat}}, \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' names(query.list)
#' alig = dtwSat(query.list[["Soybean"]], template, keep=TRUE)
#' getInternals(alig)
#' @export
setGeneric("getInternals", 
           function(object){
             object@internals
           }
)


