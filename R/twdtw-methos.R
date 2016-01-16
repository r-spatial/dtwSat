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

setMethod("show", 
          signature = signature(object="twdtw"),
          definition = function(object){
            cat("Time-Weighted DTW alignment object\n")
            cat("Number of alignments:",nrow(getAlignments(object)),"\n")
            print(head(getAlignments(object)))
            invisible(NULL)
          }
)

#' @title Summary method for twdtw-class
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Summary the alignments for each pattern in the 
#' twdtw-class object 
#' 
#' @param object An \code{\link[dtwSat]{twdtw-class}} object
#' @param ... additional arguments passed to summary
#' 
#' @docType methods
#' 
#' @return A \link[base]{data.frame} object with the the summary 
#' for each pattern 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and 
#' \code{\link[dtwSat]{twdtw}}
#'  
#' @examples
#' 
#' alig = twdtw(x=template, patterns=patterns.list)
#'        
#' show(alig)
#' summary(alig)
#' 
#' @export
setMethod("summary", 
          signature(object = "twdtw"),
          function(object, ...){
            .summary.twdtw(object, ...)
          }
)

.summary.twdtw = function(object, ...){
  res1 = do.call("rbind", lapply(object@alignments, function(pattern){
    c(N.Alig=length(pattern$distance), summary(pattern$distance))
  }))
  data.frame(res1)
}


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
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
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



###############################################################
#### GENERIC METHODS


#' @title Get pattern names from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves pattern names in the 
#' \link[dtwSat]{twdtw-class} object
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
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
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' alig = twdtw(x=template, patterns=patterns.list)
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

#' @title Get alignments from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the alignments 
#' from the object \link[dtwSat]{twdtw-class}
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
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
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun)
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



#' @title Get matching points from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the matching points 
#' for each alignment between the \code{pattern} and 
#' \code{x}
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
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
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun)
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

#' @title Get internals from twdtw object
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves cost matrix, inputs, and other 
#' internal structures from \link[dtwSat]{twdtw-class} object
#' 
#' @param object A \link[dtwSat]{twdtw-class} object
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
#'       \cr\code{x}: satellite image time series,
#'       \cr\code{N}: \code{pattern} length, and 
#'       \cr\code{M}: \code{x} length.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtw-class}}, and
#' \code{\link[dtwSat]{twdtw}}
#' 
#' @examples
#' 
#' weight.fun = logisticWeight(alpha=-0.1, beta=100, theta=0.5)
#' alig = twdtw(x=template, patterns=patterns.list, weight.fun = weight.fun, keep=TRUE)
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

