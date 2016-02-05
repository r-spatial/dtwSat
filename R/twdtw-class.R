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
#'  \item{\code{alignments}:}{A named \code{\link[base]{list}} with the same length as 
#'  the patterns list in the \code{twdtw} call. For each pattern this list stores the 
#'  results of the TWDTW analysis, that are:
#'       \cr\code{from}: a vector with the starting dates of each match,
#'       \cr\code{to}: a vector with the ending dates of each match, 
#' 	     \cr\code{distance}: a vector with TWDTW dissimilarity measure, and
#' 	     \cr\code{K}: the number of matches of the pattern.
#'  }
#' 	\item{Additional elements:}{ if \code{keep=TRUE} in the \code{twdtw} call then the list is 
#' 	extended to include internal structures used during the TWDTW computation: 
#'       \cr\code{costMatrix}: cumulative cost matrix,
#'       \cr\code{directionMatrix}: directions of steps that would be taken from each element of matrix,
#'       \cr\code{startingMatrix}: the starting points of each element of the matrix,
#'       \cr\code{stepPattern}: \code{\link[dtw]{stepPattern}} used for the 
#'       computation, see package \code{\link[dtw]{dtw}},
#'       \cr\code{N}: the length of the \code{pattern}, 
#'       \cr\code{M}: the length of the time series \code{x}, 
#'       \cr\code{timeWeight}: time weight matrix,
#'       \cr\code{localMatrix}: local cost matrix,
#'       \cr\code{patterns}: the temporal pattern, 
#'       \cr\code{x}: satellite time series,
#' 	     \cr\code{matching}: A list whose elements have the matching points for 
#' 	     each match between pattern the time series, such that:
#'          \cr--\code{index1}: a vector with matching points of the pattern, and
#'          \cr--\code{index2}: a vector with matching points of the time series.
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


