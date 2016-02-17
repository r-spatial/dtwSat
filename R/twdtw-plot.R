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

#' @title Plotting twdtw objects
#' 
#' @description Methods for plotting the results of the 
#' TWDTW analysis.
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object
#' @param type A character for the plot type: ''paths'', ''matches'', 
#' ''alignments'', ''classification'', ''cost'', ''patterns'', or ''timeseries''.
#' Default is ''alignments''. See details. 
#' @param ... additional arguments passed to plotting functions.
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
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPaths}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignments}}, 
#' \code{\link[dtwSat]{plotMatches}}, 
#' \code{\link[dtwSat]{plotClassification}},
#' \code{\link[dtwSat]{plotPatterns}}, and 
#' \code{\link[dtwSat]{plotTimeSeries}}
#' 
#' @examples 
#' 
#' log_fun = logisticWeight(alpha=-0.1, beta=100)
#' matches = twdtw(x=example_ts, patterns=patterns.list, weight.fun = log_fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' # Plot paths
#' gp1 = plot(matches, type="paths", n=1:4)
#' gp1
#' 
#' # Plot matches 
#' gp2 = plot(matches, type="matches", attr="evi")
#' gp2
#' 
#' # Plot alignments 
#' gp3 = plot(matches, type="alignments", attr=c("ndvi","evi"), threshold=4)
#' gp3
#' 
#' ## Plot classification
#' gp4 = plot(matches, type="classification", attr="evi", 
#'       from=as.Date("2009-09-01"),  to=as.Date("2013-09-01"), 
#'       by = "6 month", overlap=.5)
#' gp4
#' 
#' # Plot cost matrix
#' gp5 = plot(matches, type="cost", matrix.name="costMatrix")
#' gp5
#' 
#' # Plot cost matrix
#' gp6 = plot(matches, type="patterns")
#' gp6
#' 
#' # Plot cost matrix
#' gp7 = plot(matches, type="timeseries")
#' gp7
#' 
#' @rdname plot-method
#' 
#' @export
setMethod("plot", 
          signature(x = "twdtw"),
          function(x, type="alignments", ...){
            pt = pmatch(type,c("paths","matches","alignments","classification","cost","patterns","timeseries"))
            switch(pt,
                   plotPaths(x, ...),
                   plotMatches(x, ...),
                   plotAlignments(x, ...),
                   plotClassification(x, ...),
                   plotCostMatrix(x, ...),
                   plotPatterns(x, ...),
                   plotTimeSeries(x, ...)
            )
          }
)

