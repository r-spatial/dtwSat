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

