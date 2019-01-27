setGeneric("twdtwXtable", 
           def = function(object, ...) standardGeneric("twdtwXtable")
)

#' @title LaTeX table from accuracy metrics
#' @name twdtwXtable
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'  
#' @description Creates LaTeX table from accuracy metrics
#' 
#' @inheritParams twdtwAssessment-class
#' 
#' @param table.type Table type, 'accuracy' for User's and Producer's Accuracy, 
#' 'errormatrix' for error matrix, and 'area' for area and uncertainty. 
#' Default is 'accuracy'.
#' 
#' @param time.labels A character or numeric for the time period or NULL to 
#' include all classified periods. Default is NULL. 
#' 
#' @param category.name A character vector defining the class names. If NULL
#' the class names in the object \code{x} are used. Default is NULL.
#' 
#' @param category.type A character defining the categories type "numeric" 
#' or "letter", if NULL the class names are used. Default is NULL. 
#' 
#' @param show.prop If TRUE shows the estimated proportion of area.
#' Used with \code{table.type='accuracy'}. Default is TRUE. 
#' 
#' @param show.overall If TRUE shows the overall accuracy of the cross-validation.
#' Default is TRUE. 
#' 
#' @param rotate.col Rotate class column names in latex table. Default is FALSE. 
#' 
#' @param caption The table caption. 
#' 
#' @param digits Number of digits to show. 
#' 
#' @param conf.int Specifies the confidence level (0-1).
#' 
#' @param show.footnote Show confidence interval in the footnote. 
#' 
#' @param ... Other arguments to pass to \code{\link[xtable]{print.xtable}}.
#'
#' @seealso \code{\link[dtwSat]{twdtwAssess}} and  
#' \code{\link[dtwSat]{twdtwAssessment}}.
#'
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
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
#' # Read field samples 
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
#' twdtwXtable(twdtw_assess, table.type="errormatrix", rotate.col=TRUE,
#'   caption="Error matrix", digits=2, comment=FALSE)
#' twdtwXtable(twdtw_assess, table.type="accuracy", category.type="letter", 
#'   caption="Accuracy metrics.")
#' twdtwXtable(twdtw_assess, table.type="area", category.type="letter",
#'   digits = 0, caption="Area and uncertainty")
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
                                category.type=NULL, rotate.col=FALSE, time.labels=NULL, 
                                caption = NULL, digits = 2, show.footnote=TRUE, ...){
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
            category.colname = category.name
            if(rotate.col)
              category.colname = paste0("\\rotatebox[origin=l]{90}{",category.colname,"}")
            pt = pmatch(table.type,c("accuracy","matrix","area","errormatrix"))
            switch(pt,
                   .xtable.accuracy(x=y, category.name, category.colname, show.prop, caption, digits, show.footnote, ...),
                   .xtable.matrix(x=y, category.name, category.colname, caption, digits, ...),
                   .xtable.area(x=y, category.name, caption, digits, show.footnote, ...),
                   .xtable.matrix(x=y, category.name, category.colname, caption, digits, ...)
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
                                category.name=NULL, category.type=NULL, caption = NULL, digits = 2, show.footnote=TRUE, ...){
            y = summary(object, conf.int = conf.int)
            n = nrow(y$Users)
            if(is.null(category.name))
              category.name = rownames(y$Users)
            if(!is.null(category.type))
              category.name = switch(pmatch(category.type,c("numeric","letter")),
                                     as.character(seq(1:n)),
                                     LETTERS[1:n]
              )
            .xtable.crossvalidation(x=y, category.name, show.overall, conf.int, caption, digits, show.footnote, ...)
          }
)

.xtable.crossvalidation = function(x, category.name, show.overall, conf.int, caption, digits, show.footnote, ...){
  
  ua = sprintf(paste0("%.",digits,"f"), round(x$Users[["Accuracy"]],digits))
  ua_sd = sprintf(paste0("(%.",digits,"f)"), round(x$Users[["sd"]],digits))
  ua_ci = sprintf(paste0("[%.",digits,"f-%.",digits,"f]"), round(x$Users[["CImin"]],digits), round(x$Users[["CImax"]],digits))
  
  pa = sprintf(paste0("%.",digits,"f"), round(x$Producers[["Accuracy"]],digits))
  pa_sd = sprintf(paste0("(%.",digits,"f)"), round(x$Producers[["sd"]],digits))
  pa_ci = sprintf(paste0("[%.",digits,"f-%.",digits,"f]"), round(x$Producers[["CImin"]],digits), round(x$Producers[["CImax"]],digits))

  tbl = data.frame(ua, ua_sd, ua_ci, pa, pa_sd, pa_ci)
  table_columns = " & \\multicolumn{3}{c}{User's} & \\multicolumn{3}{c}{Producer's}"
  n = 2
  
  if(show.overall){
    oa = sprintf(paste0("%.",digits,"f"), round(x$Overall[["Accuracy"]],digits))
    oa_sd = sprintf(paste0("(%.",digits,"f)"), round(x$Overall[["sd"]],digits))
    oa_ci = sprintf(paste0("[%.",digits,"f-%.",digits,"f]"), round(x$Overall[["CImin"]],digits), round(x$Overall[["CImax"]],digits))

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
                              "\\multicolumn{1}{c}{Class}", paste(rep(" & \\multicolumn{1}{c}{$\\mu$} & \\multicolumn{1}{c}{$\\sigma$} & \\multicolumn{1}{c}{ci*}", n),collapse = ""),"\\\\\n"),
                       paste("\\hline \n", ifelse(show.footnote, paste0("\\multicolumn{",ncol(tbl)+1,"}{l}{* ",conf.int*100,"\\% confidence interval.}\n"), ""), sep = ""))
  
  rownames(tbl) = category.name
  
  tbl = xtable(tbl, caption)
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0), sanitize.text.function = function(x) x, ...)
  
}


.xtable.accuracy = function(x, category.name, category.colname, show.prop, caption, digits, show.footnote, ...){

  ua = sprintf(paste0("%.",digits,"f$\\pm$%.",digits,"f"), round(x$UsersAccuracy[,"Accuracy"],digits), round(x$UsersAccuracy[,"ci"], digits))
  pa = sprintf(paste0("%.",digits,"f$\\pm$%.",digits,"f"), round(x$ProducersAccuracy[,"Accuracy"],digits), round(x$ProducersAccuracy[,"ci"], digits))
  oa = c(sprintf(paste0("%.",digits,"f$\\pm$%.",digits,"f"), round(x$OverallAccuracy["Accuracy"],digits), round(x$OverallAccuracy["ci"], digits)), rep("", length(pa)-1))
  
  prop = data.frame(ua, pa, oa)
  names(prop) = c("User's*", "Producers's*", "Overall*")
  
  if(show.prop){
    prop = data.frame(`User's*` = "", `Producers's*` = "", `Overall*` = "")
    prop = as.data.frame.matrix(x$ProportionMatrix)
    prop = data.frame(apply(prop[,!names(prop)%in%c("Area","w")], 1, FUN = sprintf, fmt=paste0(paste0("%.",digits,"f"))), stringsAsFactors = FALSE)
    rownames(prop) = names(prop)    
    prop$`User's*` = ""
    prop$`Producers's*` = ""
    prop$`Overall*` = ""
    names(prop)[1:length(category.name)] = category.colname
    prop$`User's*`[1:length(ua)] = ua
    prop$`Producers's*`[1:length(pa)] = pa
    prop$`Overall*`[1:length(oa)] = oa
  }

  rownames(prop)[1:length(category.name)] = category.name
  tbl = xtable(prop, caption)
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$pos[[2]] = c(nrow(tbl))
  if(show.prop){
    comment$command  = c(paste0("&\\multicolumn{",length(category.name),"}{c}{Reference class}&&&&\\\\\n", 
                                paste(c("\\multicolumn{1}{c}{Map class}",names(tbl)), collapse = " & "),"\\\\\n"),
                         paste("\\hline \n", ifelse(show.footnote, paste0("\\multicolumn{",ncol(tbl),"}{l}{* ",x$conf.int*100,"\\% confidence interval.}\n"), ""), sep = ""))
  } else {
    comment$command  = c(paste0(paste(c("\\multicolumn{1}{c}{Class}",names(tbl)), collapse = " & "),"\\\\\n"),
                         paste("\\hline \n", ifelse(show.footnote, paste0("\\multicolumn{",ncol(tbl),"}{l}{* ",x$conf.int*100,"\\% confidence interval.}\n"), ""), sep = ""))
  }

  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
                       hline.after = c(-1, 0), sanitize.text.function = function(x) x, ...)
}

.xtable.matrix = function(x, category.name, category.colname, caption, digits, ...){
  m = as.data.frame.matrix(x$ErrorMatrix)
  # names(m)[ncol(m)] = "Estimation weight"
  names(m)[1:length(category.name)] = category.colname
  rownames(m)[1:length(category.name)] = category.name
  
  tbl = xtable(m, caption, digits = c(rep(0, ncol(m)-1), digits, 3))

  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$command  = c(paste0("&\\multicolumn{",length(category.name),"}{c}{Reference class}&&\\\\\n", 
                              paste(c("\\multicolumn{1}{c}{Map class}",names(tbl)), collapse = " & "),"\\\\\n"))
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0, nrow(tbl)), sanitize.text.function = function(x) x, ...)
  
}

.xtable.area = function(x, category.name, caption, digits, show.footnote, ...){
  
  a = x$AreaUncertainty
  a = data.frame(a)

  mp = sprintf(paste0("%.",digits,"f"), round(unlist(a$Mapped),digits))
  ad = sprintf(paste0("%.",digits,"f"), round(unlist(a$Adjusted),digits))
  ci = sprintf(paste0("$\\pm$%.",digits,"f"), round(unlist(a$ci),digits))
  
  tbl = data.frame(mp, ad, ci)
  rownames(tbl) = category.name
  names(tbl) = c("Mapped area", "Adjusted area", "Margin of error*")
  tbl = xtable(tbl, caption)
  
  comment          = list()
  comment$pos      = list()
  comment$pos[[1]] = c(0)
  comment$pos[[2]] = c(nrow(tbl))
  comment$command  = c(paste0(paste(c("\\multicolumn{1}{c}{Class}",names(tbl)), collapse = " & "), "\\\\\n"),
                       paste("\\hline \n", ifelse(show.footnote, paste0("\\multicolumn{",ncol(tbl),"}{l}{* ",x$conf.int*100,"\\% confidence interval.}\n"), ""), sep = ""))
  
  print.xtable(tbl, add.to.row = comment, include.rownames=TRUE, include.colnames = FALSE,
               hline.after = c(-1, 0), sanitize.text.function = function(x) x, ...)
  
}



