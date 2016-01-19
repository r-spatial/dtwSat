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
#   R Package dtwSat - 2016-16-01                             #
#                                                             #
###############################################################

#' @title Plotting twdtw objects
#' 
#' @description Methods for plotting the results of the 
#' Time-Weighted DTW analysis
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @param x A \code{\link[dtwSat]{twdtw-class}} object
#' @param type A character for the plot type, ''path'', ''match'', 
#' ''alignment'', ''group'', ''cost''. Default is ''path''
#' @param ... additional arguments passed to plotting functions
#' \code{\link[dtwSat]{twdtw-class}}, \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}}, \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @return A \link[ggplot2]{ggplot} object
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, 
#' \code{\link[dtwSat]{twdtw}}, 
#' \code{\link[dtwSat]{plotPath}}, 
#' \code{\link[dtwSat]{plotCostMatrix}},
#' \code{\link[dtwSat]{plotAlignment}},
#' \code{\link[dtwSat]{plotMatch}}, and
#' \code{\link[dtwSat]{plotGroup}}
#' 
#' @examples 
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, 
#'         normalize.patterns=TRUE, patterns.length=23, keep=TRUE)
#' 
#' # Plot paths
#' gp1 = plot(alig, type="path", n.alignments=1:4)
#' gp1
#' 
#' # Plot matches 
#' gp2 = plot(alig, type="match", attr="evi")
#' gp2
#' 
#' # Plot alignments 
#' gp3 = plot(alig, type="alignment", attr=c("ndvi","evi"), threshold=4)
#' gp3
#' 
#' ## Plot classification
#' gp4 = plot(alig, type="group", attr="evi", from=as.Date("2009-09-01"),  
#'            to=as.Date("2013-09-01"), by = "6 month", overlap=.3)
#' gp4
#' 
#' # Plot cost matrix
#' gp5 = plot(alig, type="cost", matrix.name="costMatrix")
#' gp5
#' 
#' # Plot cost matrix
#' gp6 = plot(alig, type="pattern")
#' gp6
#' 
#' @export
setMethod("plot", 
          signature(x = "twdtw"),
          function(x, type="path", ...){
            if(!is(x,"twdtw"))
              stop("x is not a twdtw object.")
            if(length(getInternals(x))==0)
              stop("plot method requires twdtw internals (set keep.internals=TRUE on dtw() call)")
            pt = pmatch(type,c("path","match","alignment","group","cost","pattern"))
            switch(pt,
                   plotPath(x, ...),
                   plotMatch(x, ...),
                   plotAlignment(x, ...),
                   plotGroup(x, ...),
                   plotCostMatrix(x, ...),
                   plotPatterns(x, ...)
            )
          } 
)

