setGeneric("twdtwXtable", 
           def = function(object, ...) standardGeneric("twdtwXtable")
)

#' @title Latex table from accuracy metrics
#' @name twdtwXtable
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'  
#' @description Creates Latex table from accuracy metrics
#' 
#' @inheritParams twdtwAssessment-class
#' 
#' @param table.type table type, 'accuracy' for User's and Producer's Accuracy, 
#' 'errormatrix' for error matrix, and 'area' for area and uncertainty. 
#' Default is 'accuracy'.
#' 
#' @param time.labels a character or numeric for the time period or NULL to 
#' include all classified periods. Default is NULL. 
#' 
#' @param category.name a character vector defining the class names. If NULL
#' then use the classe names in the object \code{x}. Default is NULL.
#' 
#' @param category.type a character defining the categories type "numeric" 
#' or "letter", if NULL then use the class names. Default is NULL. 
#' 
#' @param show.prop if TRUE shows the estimated proportion of area.
#' Used with \code{table.type='accuracy'}. Default is TRUE. 
#' 
#' @param show.overall if TRUE shows the overall accuracy of the cross-validation.
#' Default is TRUE. 
#' 
#' @param conf.int specifies the confidence level (0-1).
#' 
#' @param ... other arguments to pass to \code{\link[xtable]{xtable}}.
#'
#' @seealso \code{\link[dtwSat]{twdtwAssess}} and  
#' \code{\link[dtwSat]{twdtwAssessment}}.
#'
#' @examples 
#' \dontrun{
#' 
#' # Create raster time series
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' 
#' # Read fiels samples 
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", 
#'                 package="dtwSat"), what = "character")
#' 
#' # Split samples for training (10%) and validation (90%) using stratified sampling 
#' library(caret) 
#' set.seed(1)
#' I = unlist(createDataPartition(field_samples$label, p = 0.1))
#' training_samples = field_samples[I,]
#' validation_samples = field_samples[-I,]
#' 
#' # Create temporal patterns 
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))
#' 
#' # Run TWDTW analysis for raster time series 
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' r_twdtw = twdtwApply(x=rts, y=temporal_patterns, weight.fun=log_fun, format="GTiff", 
#'                      overwrite=TRUE)
#'                      
#' # Classify raster based on the TWDTW analysis 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' plot(r_lucc)
#' 
#' # Assess classification 
#' twdtw_assess = twdtwAssess(object = r_lucc, y = validation_samples, 
#'                            proj4string = proj_str, conf.int=.95) 
#' twdtw_assess
#' 
#' # Create latex tables 
#' twdtwXtable(twdtw_assess, table.type="errormatrix", category.type="letter")
#' twdtwXtable(twdtw_assess, table.type="accuracy", category.type="letter")
#' twdtwXtable(twdtw_assess, table.type="area", category.type="letter")
#' 
#' }
NULL

#' @aliases twdtwXtable-twdtwAssessment
#' @inheritParams twdtwXtable
#' @rdname twdtwXtable 
#' @export
setMethod("twdtwXtable", 
          signature = signature(object = "twdtwAssessment"),
          definition = function(object, table.type="accuracy", show.prop=TRUE, category.name=NULL,
                                category.type=NULL, time.labels=NULL, ...){
            y = object@accuracySummary
            if(!is.null(time.labels))
              y = object@accuracyByPeriod[[time.labels]]
            if(is.null(y))
              stop("time.labels out of bounds", call. = TRUE)
            n = nrow(object@accuracySummary$ProportionMatrix) - 1
            if(is.null(category.name))
              category.name = rownames(object@accuracySummary$ProportionMatrix)[-(n+1)]
            if(!is.null(category.type))
              category.name = switch(pmatch(category.type,c("numeric","letter")),
                   as.character(seq(1:n)),
                   LETTERS[1:n]
              )
            pt = pmatch(table.type,c("accuracy","matrix","area","errormatrix"))
            switch(pt,
                   .xtable.accuracy(x=y, category.name, show.prop, ...),
                   .xtable.matrix(x=y, category.name, ...),
                   .xtable.area(x=y, category.name, ...),
                   .xtable.matrix(x=y, category.name, ...)
            )
          }
)

#' @aliases twdtwXtable-twdtwCrossValidation
#' @inheritParams twdtwXtable
#' @rdname twdtwXtable 
#' @export
setMethod("twdtwXtable", 
          signature = signature(object = "twdtwCrossValidation"),
          definition = function(object, conf.int=.95, show.overall=TRUE, 
                                category.name=NULL, category.type=NULL, ...){
            y = summary(object, conf.int = conf.int)
            n = nrow(y$Users)
            if(is.null(category.name))
              category.name = rownames(y$Users)
            if(!is.null(category.type))
              category.name = switch(pmatch(category.type,c("numeric","letter")),
                                     as.character(seq(1:n)),
                                     LETTERS[1:n]
              )
            .xtable.crossvalidation(x=y, category.name, show.overall, conf.int, ...)
          }
)

.xtable.crossvalidation = function(x, category.name, show.overall, conf.int, ...){
  
  ua = sprintf("%.2f", round(x$Users[["Accuracy"]],2))
  ua_sd = sprintf("(%.2f)", round(x$Users[["sd"]],2))
  ua_ci = sprintf("[%.2f-%.2f]", round(x$Users[["CImin"]],2), round(x$Users[["CImax"]],2))
  
  pa = sprintf("%.2f", round(x$Producers[["Accuracy"]],2))
  pa_sd = sprintf("(%.2f)", round(x$Producers[["sd"]],2))
  pa_ci = sprintf("[%.2f-%.2f]", round(x$Producers[["CImin"]],2), round(x$Producers[["CImax"]],2))

  tbl = data.frame(ua, ua_sd, ua_ci, pa, pa_sd, pa_ci)
  table_columns = " & \\multicolumn{3}{c}{User's} & \\multicolumn{3}{c}{Producer's}"
  n = 2
  
  if(show.overall){
    oa = sprintf("%.2f", round(x$Overall[["Accuracy"]],2))
    oa_sd = sprintf("(%.2f)", round(x$Overall[["sd"]],2))
    oa_ci = sprintf("[%.2f-%.2f]", round(x$Overall[["CImin"]],2), round(x$Overall[["CImax"]],2))

    tbl$oa = ""
    tbl$oa_sd = ""
    tbl$oa_ci = ""
    
    tbl$oa[1] = oa
    tbl$oa_sd[1] = oa_sd
    tbl$oa_ci[1] = oa_ci 
    
    table_columns = paste0(table_columns, " & \\multicolumn{3}{c}{Overall}")
    n = 3
  }
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$pos[[2]] = c(nrow(tbl))
  comment$command  = c(paste0(table_columns, "\\\\\n", 
                              "Class", paste(rep(" & $\\mu$ & $\\sigma$ & ci*", n),collapse = ""),"\\\\\n"),
                       paste("\\hline \n", "\\multicolumn{",ncol(tbl)+1,"}{l}{* ",conf.int*100,"\\% confidence interval.}\n", sep = ""))
  
  rownames(tbl) = category.name
  
  tbl = xtable(tbl, ...)
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0), sanitize.text.function = function(x) x)
  
}


.xtable.accuracy = function(x, category.name, show.prop, ...){
  
  prop = x$ProportionMatrix
  prop = data.frame(apply(prop[,!names(prop)%in%c("Area","w")], 1, FUN = sprintf, fmt="%.2f"), stringsAsFactors = FALSE)
  rownames(prop) = names(prop)
  prop$`User's*` = ""
  prop$`Producers's*` = ""
  prop$`Overall*` = ""

  ua = sprintf("%.2f$\\pm$%.2f", round(x$UsersAccuracy[,"Accuracy"],2), round(x$UsersAccuracy[,"ci"], 2))
  pa = sprintf("%.2f$\\pm$%.2f", round(x$ProducersAccuracy[,"Accuracy"],2), round(x$ProducersAccuracy[,"ci"], 2))
  oa = sprintf("%.2f$\\pm$%.2f", round(x$OverallAccuracy["Accuracy"],2), round(x$OverallAccuracy["ci"], 2))
  
  prop$`User's*`[1:length(ua)] = ua
  prop$`Producers's*`[1:length(pa)] = pa
  prop$`Overall*`[1:length(oa)] = oa
  names(prop)[1:length(category.name)] = category.name
  rownames(prop)[1:length(category.name)] = category.name
  tbl = xtable(prop, ...)
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$pos[[2]] = c(nrow(tbl))
  comment$command  = c(paste0("&\\multicolumn{",ncol(tbl)-1,"}{c}{Reference class}\\\\\n", 
                              paste(c("Map class",names(tbl)), collapse = " & "),"\\\\\n"),
                       paste("\\hline \n", "\\multicolumn{",ncol(tbl),"}{l}{* ",x$conf.int*100,"\\% confidence interval.}\n", sep = ""))
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
                       hline.after = c(-1, 0), sanitize.text.function = function(x) x)
}

.xtable.matrix = function(x, category.name, ...){
  m = x$ErrorMatrix
  names(m)[ncol(m)] = "Estimation weight"
  names(m)[1:length(category.name)] = category.name
  rownames(m)[1:length(category.name)] = category.name
  
  tbl = xtable(m, digits = c(rep(0, ncol(m)-1), 2, 2), ...)
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$command  = c(paste0("&\\multicolumn{",ncol(tbl)-1,"}{c}{Reference class}\\\\\n", 
                              paste(c("Map class",names(tbl)), collapse = " & "),"\\\\\n"))
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0), sanitize.text.function = function(x) x)
  
}

.xtable.area = function(x, category.name, ...){
  
  a = x$AreaUncertainty
  a = data.frame(a)
  
  mp = sprintf("%.2f", round(unlist(a$Mapped),2))
  ad = sprintf("%.2f", round(unlist(a$Adjusted),2))
  ci = sprintf("$\\pm$%.2f", round(unlist(a$ci),2))
  
  tbl = data.frame(mp, ad, ci)
  rownames(tbl) = category.name
  names(tbl) = c("Mapped area", "Adjusted area", "Margin of error*")
  tbl = xtable(tbl, ...)
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$pos[[2]] = c(nrow(tbl))
  comment$command  = c(paste0(paste(c("Class",names(tbl)), collapse = " & "), "\\\\\n"),
                       paste("\\hline \n", "\\multicolumn{",ncol(tbl),"}{l}{* ",x$conf.int*100,"\\% confidence interval.}\n", sep = ""))
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0), sanitize.text.function = function(x) x)
  
}



