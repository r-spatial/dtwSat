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
#' ''alignments'', ''classification'', ''cost'', ''patterns'', or ''timeseries''.
#' 
#' @param ... additional arguments to pass to plotting functions.
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotMatches}}, 
#' \code{\link[dtwSat]{plotClassification}},
#' \code{\link[dtwSat]{plotPatterns}}, or 
#' \code{\link[dtwSat]{plotTimeSeries}}.
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
#' @examples 
#' 
#' # Plot patterns 
#' patt = twdtwTimeSeries(patterns.list)
#' gp1 = plot(patt, type="patterns")
#' gp1
#' 
#' ts = twdtwTimeSeries(example_ts)
#' gp2 = plot(ts, type="timeseries")
#' gp2
#' 
#' @export
NULL


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


