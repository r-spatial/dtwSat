#' @aliases xtable
#' @inheritParams twdtwAssessment-class
#' @rdname twdtwAssessment-class
#' 
#' @param x an object of class \code{\link[dtwSat]{twdtwAssessment}}.
#' @param type table type, 'accuracy' for User's and Producer's Accuracy, 
#' 'matrix' for error matrix, and 'area' for area and uncertainty. 
#' Default is 'accuracy'.
#' @param time.labels a character or numeric for the time period or NULL to 
#' include all classified periods. Default is NULL. 
#' @param ... other arguments to pass to \code{\link[xtable]{xtable}}.
#' 
#' @export
setMethod("xtable", 
          signature = signature(x = "twdtwAssessment"),
          definition = function(x, type="accuracy", time.labels=NULL, ...){
            xtable(x@accuracySummary$ErrorMatrix, ...)
            pt = pmatch(type,c("accuracy","matrix","area"))
            switch(pt,
                   .xtable.accuracy(x=x@accuracySummary, ...),
                   .xtable.matrix(x=x@accuracySummary$ErrorMatrix, ...),
                   .xtable.area(x=x@accuracySummary$AreaUncertainty, ...)
            )
          }
)

.xtable.accuracy = function(x, ...){

  prop = x$ProportionMatrix
  prop = data.frame(apply(prop[,!names(prop)%in%c("Area","w")], 1, FUN = sprintf, fmt="%.2f"), stringsAsFactors = FALSE)
  prop$`User's` = ""
  prop$`Producers's` = ""
  prop$`Overall` = ""

  ua = sprintf("%.2f$\\pm$%.2f", round(x$UsersAccuracy[,"Accuracy"],2), round(x$UsersAccuracy[,"ci"], 2))
  pa = sprintf("%.2f$\\pm$%.2f", round(x$ProducersAccuracy[,"Accuracy"],2), round(x$ProducersAccuracy[,"ci"], 2))
  oa = sprintf("%.2f$\\pm$%.2f", round(x$OverallAccuracy["Accuracy"],2), round(x$OverallAccuracy["ci"], 2))
  
  prop$`User's`[1:length(ua)] = ua
  prop$`Producers's`[1:length(pa)] = pa
  prop$`Overall`[1:length(oa)] = oa
  tbl = xtable(prop)
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(nrow(tbl))
  comment$command  = c(paste("\\hline \n",  "Test foot note.  \n", sep = ""))
  
  print.xtable(tbl, add.to.row = comment, sanitize.text.function = function(x) x)
}

.xtable.matrix = function(x, ...){
  tbl = xtable(x, digits = c(rep(0, ncol(x)-1), 2, 2), ...)
  print.xtable(tbl, sanitize.text.function = function(x) x)
}

.xtable.area = function(x, ...){
  
  x=twdtw_assess@accuracySummary$AreaUncertainty
  x = data.frame(x)
  
  mp = sprintf("%.2f", round(unlist(x$Mapped),2))
  ad = sprintf("%.2f", round(unlist(x$Adjusted),2))
  ci = sprintf("$\\pm$%.2f", round(unlist(x$ci),2))
  
  tbl = data.frame(mp, ad, ci)
  names(tbl) = c("Mapped", "Adjusted", "Margin of error (95\\% CI)")
  
  tbl = xtable(tbl)
  
  print.xtable(tbl, sanitize.text.function = function(x) x)
}



