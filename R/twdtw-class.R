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


#' @title twdtw-class
#' @name twdtw-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'
#' @description Class for Time-Weighted Dynamic Time Warping matches.
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
#'       \cr\code{M}: \code{x} length, 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{patterns}: temporal pattern, 
#'       \cr\code{x}: satellite image time series,
#'       }
#' 	\item{\code{matching}:}{An object of class \code{\link[base]{list}} whose 
#'   elements have the matching points for each alignment between the 
#'   \code{pattern} and the \code{x}, such that
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
NULL
twdtw = setClass(
  Class = "twdtw",
  slots = c(
    call = "call",
    alignments = "list"
  ),
  validity = function(object){
    if(!is(object@call, "call")){
      stop("[twdtw: validation] Invalid call object, class different from call.")
    }else{}
    if(!is(object@alignments, "list")){
      stop("[twdtw: validation] Invalid TWDTW alignments, class different from data.frema.")
    }else{}
    return(TRUE)
  }
)

setMethod("initialize",
  signature = "twdtw",
  definition = 
    function(.Object, call, alignments){
      .Object@call = new("call")
      .Object@alignments = list(list(
                                from       = numeric(0), 
                                to         = numeric(0), 
                                distance   = numeric(0),
                                K          = numeric(0),
                                pattern      = numeric(0),
                                x = numeric(0),
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


