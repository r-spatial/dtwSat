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
#'       \cr\code{pattern}: pattern identification,
#'       \cr\code{from}: starting date,
#'       \cr\code{to}: ending date, 
#' 	     \cr\code{distance}: TWDTW distance, and
#' 	     \cr\code{K}: the number of alignments of the pattern.
#'       }
#' 	\item{\code{internals}:}{An object of class \code{\link[base]{list}} whose 
#'   elements have the internal structures used by \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken in the alignments,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}}
#'       \cr\code{N}: \code{pattern} length, 
#'       \cr\code{M}: \code{timeseries} length, 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{patterns}: temporal pattern, 
#'       \cr\code{timeseries}: satellite image time series,
#'       }
#' 	\item{\code{matching}:}{An object of class \code{\link[base]{list}} whose 
#'   elements have the matching points for each alignment between the 
#'   \code{pattern} and the \code{timeseries}, such that
#'       \cr\code{index1}: matching points of the temporal pattern, and
#'       \cr\code{index2}: matching points of the time series.
#'       }
#' }
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{getPatternNames}},
#' \code{\link[dtwSat]{getAlignments}},
#' \code{\link[dtwSat]{getMatches}}, and 
#' \code{\link[dtwSat]{getInternals}}
#' 
NULL
dtwSat = setClass(
  Class = "dtwSat",
  slots = c(
    call = "call",
    alignments = "list"
  ),
  validity = function(object){
    if(!is(object@call, "call")){
      stop("[dtwSat: validation] Invalid call object, class different from call.")
    }else{}
    if(!is(object@alignments, "list")){
      stop("[dtwSat: validation] Invalid TWDTW alignments, class different from data.frema.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "dtwSat",
  definition = 
    function(.Object, call, alignments){
      .Object@call = new("call")
      .Object@alignments = list(list(
                                from       = numeric(0), 
                                to         = numeric(0), 
                                distance   = numeric(0),
                                K          = numeric(0),
                                pattern      = numeric(0),
                                timeseries = numeric(0),
                                matching   = list(
                                                index1 = numeric(0), 
                                                index2 = numeric(0)
                                                ),
                                internals  = list(0)
                                ))
      if(!missing(call))
        .Object@call = call
      if(!missing(alignments))
        .Object@alignments = alignments
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

setMethod("summary", 
          signature(object = "dtwSat"),
          function(object, ...){
            res1 = do.call("rbind", lapply(object@alignments, function(pattern){
              c(N.Alig=length(pattern$distance), summary(pattern$distance))
            }))
            class(res1) = c("summaryDefault", "table", oldClass(object))
            res1
          }
)


###############################################################
#### GENERIC METHODS


#' @title Get pattern names from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves pattern names in the 
#' \link[dtwSat]{dtwSat-class} object
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @param p.names A \link[base]{character} or \link[base]{numeric}
#' vector with the patterns identification. If not declared the function 
#' retrieves the names for all patterns 
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{character}}
#' or \code{\link[base]{numeric}} vector 
#' 
#' @seealso 
#' \code{\link[dtwSat]{dtwSat-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' alig = twdtw(patterns=patterns.list, timeseries=template)
#' 
#' getPatternNames(alig)
#' 
#' getPatternNames(alig, p.names=c(1,3))
#' 
#' getPatternNames(alig, p.names="Maize")
#' 
#' @export
setGeneric("getPatternNames", 
           function(object, p.names){
             p.names = .getPatternNames(object, p.names)
             if(any(is.na(p.names)))
               warning("the patterns identification is invalid", call. = FALSE)
             p.names
           }
)

.getPatternNames = function(object, p.names){
  if(missing(p.names)) p.names = seq_along(object@alignments)
  all_names = names(object@alignments)
  names(all_names) = all_names
  if(is.null(all_names)) all_names = seq_along(object@alignments)
  all_names[p.names]
}

#' @title Get alignments from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the alignments 
#' from the object \link[dtwSat]{dtwSat-class}
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \link[base]{data.frame} with the following columns:  
#'       \cr\code{pattern}: the pattern identification,
#'       \cr\code{from}: starting date,
#'       \cr\code{to}: ending date, and
#' 	     \cr\code{distance}: TWDTW distances.
#' 
#' @seealso 
#' \code{\link[dtwSat]{dtwSat-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(patterns=patterns.list, timeseries=template, weight.fun = weight.fun)
#' 
#' getAlignments(alig)
#' 
#' getAlignments(alig, p.names="Soybean")
#' 
#' getAlignments(alig, p.names=c(2,3))
#' 
#' @export
setGeneric("getAlignments", 
           function(object, ...) {
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getAlignments(object, p.names)
           }
)

.getAlignments = function(object, p.names){
  k = 1
  names(p.names) = NULL
  do.call("rbind", lapply(p.names, function(p){
    name = p
    x = object@alignments[[p]]
    if(length(x$distance)<1)
      name = numeric(0)
    r.names = paste0(k:(k+x$K-1))
    k <<- k + x$K
    data.frame(pattern  = name,
               from     = x$from,
               to       = x$to,
               distance = x$distance,
               stringsAsFactors = FALSE,
               row.names = r.names
               )
  }))
}



#' @title Get matching points from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the matching points 
#' for each alignment between the \code{pattern} and the 
#' \code{timeseries}
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} whose 
#'  elements have the matching points for each alignment between 
#'  the temporal pattern and the time series. 
#'  Each element has two vectors: 
#'       \cr\code{index1}: matching points of the pattern, and
#'       \cr\code{index2}: matching points of the time series.
#' 
#' @seealso 
#' \code{\link[dtwSat]{dtwSat-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(patterns=patterns.list, timeseries=template, weight.fun = weight.fun)
#' 
#' getMatches(alig)
#' 
#' getMatches(alig, p.names="Maize")
#' 
#' getMatches(alig, p.names=1)
#' 
#' @export
setGeneric("getMatches", 
           function(object, ...){
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getMatches(object, p.names)
           }
)

.getMatches = function(object, p.names){
  lapply(p.names, function(p) object@alignments[[p]]$matching)
}

#' @title Get internals from dtwSat object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves cost matrix, inputs, and other 
#' internal structures from \link[dtwSat]{dtwSat-class} object
#' 
#' @param object A \link[dtwSat]{dtwSat-class} object
#' @param ... additional arguments passed to \code{\link[dtwSat]{getPatternNames}}
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} whose 
#'   elements have the internal structures used used in \code{\link[dtwSat]{twdtw}}. 
#'   The elements are: 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken in the alignments,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}}
#'       \cr\code{pattern}: pattern time series, 
#'       \cr\code{timeseries}: satellite image time series,
#'       \cr\code{N}: \code{pattern} length, and 
#'       \cr\code{M}: \code{timeseries} length.
#' 
#' @seealso 
#' \code{\link[dtwSat]{dtwSat-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(patterns=patterns.list, timeseries=template, weight.fun = weight.fun, keep=TRUE)
#' 
#' a = getInternals(alig)
#' names(a) 
#' 
#' a = getInternals(alig, p.names="Maize")
#' names(a) 
#' 
#' a = getInternals(alig, p.names=c(1,2))
#' names(a) 
#' 
#' 
#' @export
setGeneric("getInternals", 
           function(object, ...){
             p.names = getPatternNames(object, ...)
             if(any(is.na(p.names)))
               stop("the patterns identification is invalid")
             .getInternals(object, p.names)
           }
)

.getInternals = function(object, p.names) {
  lapply(p.names, function(p) object@alignments[[p]]$internals)
}








