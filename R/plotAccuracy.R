#' @title Plotting accuracy assessment 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting accuracy assessment results.
#' 
#' @param x An object of class \code{\link[dtwSat]{twdtwAssessment}} or 
#' \code{\link[dtwSat]{twdtwCrossValidation}}.
#' 
#' @param time.labels a character or numeric for the time periods or NULL to 
#' aggregate all classified periods in the same plot. Default is NULL. Used 
#' if \code{x} is \code{\link[dtwSat]{twdtwAssessment}}.
#' 
#' @param perc if TRUE shows the results in percent of area. Otherwise shows the 
#' area in the map units or km2 for no project raster. Default is TRUE.
#' 
#' @param category.name a character vector defining the class names. If NULL
#' the class names in the object \code{x} are used. Default is NULL.
#' 
#' @param category.type a character defining the categories type "numeric" 
#' or "letter", if NULL the class names are used. Default is NULL. 
#' 
#' @param conf.int confidence level (0-1) for interval estimation of the population mean.
#' For details see \code{\link[Hmisc]{smean.cl.normal}}. Used if \code{x} is 
#' \code{\link[dtwSat]{twdtwCrossValidation}}.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwAssessment}} and \code{\link[dtwSat]{twdtwAssess}}
#'  
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples
#' \dontrun{
#' 
#' # See ?twdtwAssess and ?twdtwCrosValidate
#' 
#' plotAccuracy(x)
#' 
#' plotAccuracy(x, category.type="letter")
#' 
#' }
#' 
#' @export
plotAccuracy = function(x, perc=TRUE, conf.int=.95, time.labels=NULL, 
                        category.name=NULL, category.type=NULL){
  
  if(class(x)=="twdtwCrossValidation"){
    gp = .plotCrossValidation(x, conf.int, perc, category.name, category.type)
  } else {
    if(class(x)=="twdtwAssessment"){
      gp = .plotAssessmentAccuracy(x, perc, time.labels, category.name, category.type)
    } else {
      stop("Class of x is not twdtwAssessment or twdtwCrossValidation")
    }
  } 
  
  gp 
  
}

.plotAssessmentAccuracy = function(x, perc=TRUE, time.labels=NULL, 
                            category.name=NULL, category.type=NULL){
  
  if(is.null(category.name)){
    category.name = rownames(x@accuracySummary$ProportionMatrix)
    category.name = category.name[-length(category.name)]
  }
  if(!is.null(category.type))
    category.name = switch(pmatch(category.type,c("numeric","letter")),
                           as.character(seq(1:length(category.name))),
                           LETTERS[1:length(category.name)]
    )
  
  y = list(`Accumulated` = x@accuracySummary)
  if(!is.null(time.labels))
    y = x@accuracyByPeriod[time.labels]
  if(is.null(y))
    stop("time.labels out of bounds", call. = TRUE)

  time.labels = names(y)
  if(is.null(time.labels))
    time.labels = seq_along(y)
  
  df = do.call("rbind", lapply(time.labels, function(i){
    # User's 
    df           = data.frame(y[[i]]$UsersAccuracy)
    df_u         = data.frame(value    = df$Accuracy,
                              variable = category.name,
                              Legend   = "User's",
                              ci       = df$ci,
                              Period  = i)
    # Producer's
    df           = data.frame(y[[i]]$ProducersAccuracy)
    df_p         = data.frame(value    = df$Accuracy,
                              variable = category.name,
                              Legend   = "Producer's",
                              ci       = df$ci,
                              Period  = i)
    
    df           = rbind(df_u, df_p)
    df
  }))
  
  limits = aes_string(ymax = "value + ci", ymin = "value - ci")
  dodge = position_dodge(width=0.9)
  
  gp = ggplot(df, aes_string(fill="Legend", y="value", x="variable")) + 
    facet_wrap(~Period, scales = "free") + 
    geom_bar(position="dodge", stat="identity", na.rm=TRUE) +
    geom_errorbar(limits, position=dodge, width=0.25, na.rm=TRUE) + 
    xlab("") + 
    ylab("Accuracy")
  
  if(perc)
    gp = gp + scale_y_continuous(labels = percent)
  
  gp
  
}

.plotCrossValidation = function(x, conf.int, perc, category.name, category.type){
  
  if(is.null(category.name)){
    category.name = rownames(x@accuracy$Resample1$ErrorMatrix)
    category.name = category.name
  }
  if(!is.null(category.type))
    category.name = switch(pmatch(category.type,c("numeric","letter")),
                           as.character(seq(1:length(category.name))),
                           LETTERS[1:length(category.name)]
    )
  
  UA = do.call("rbind", lapply(x@accuracy, function(x) data.frame(label="UA", rbind(x$UsersAccuracy))))
  names(UA)[-1] = category.name
  PA = do.call("rbind", lapply(x@accuracy, function(x) data.frame(label="PA", rbind(x$ProducersAccuracy))))
  names(PA)[-1] = category.name
  df = melt(rbind(UA,PA), id="label")
  df$label = factor(df$label, levels = c("UA", "PA"), 
                    labels = c("User's Accuracy", "Producer's Accuracy"))
  df$variable = factor(df$variable, levels = levels(df$variable), 
                       labels = gsub("[.]","-",levels(df$variable)))
  
  gp = ggplot(df, aes_string(x="variable", y="value")) +
    stat_summary(fun.data="mean_cl_boot", fun.args=list(conf.int = conf.int), 
                 width=0.5, geom="crossbar", size=0.1, fill = "gray") + 
    geom_point(size=0.2) +  
    facet_grid(. ~ label) + 
    xlab("") + 
    ylab("Accuracy") + 
    coord_flip()
  
  if(perc){
    gp = gp + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2), labels = percent)
  } else {
    gp = gp + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2)) 
  }
    
  
  gp
}

