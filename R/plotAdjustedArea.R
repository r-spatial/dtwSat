#' @title Plotting area and uncertainty 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting area and uncertainty.
#' 
#' @inheritParams plotAccuracy
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
#' # See ?twdtwAssess
#' 
#' plotAdjustedArea(twdtw_assess)
#' 
#' plotAdjustedArea(twdtw_assess, category.type="letter")
#' 
#' }
#' 
#' @export
plotAdjustedArea = function(x, perc=TRUE, time.labels=NULL, 
                            category.name=NULL, category.type=NULL){
  
  y = list(`Accumulated area` = x@accuracySummary)
  if(!is.null(time.labels))
    y = x@accuracyByPeriod[time.labels]
  if(is.null(y))
    stop("time.labels out of bounds", call. = TRUE)

  if(is.null(category.name)){
    category.name = rownames(x@accuracySummary$ProportionMatrix)
    category.name = category.name[-length(category.name)]
  }
  if(!is.null(category.type))
    category.name = switch(pmatch(category.type,c("numeric","letter")),
                           as.character(seq(1:length(category.name))),
                           LETTERS[1:length(category.name)]
    )
  
  time.labels = names(y)
  if(is.null(time.labels))
    time.labels = seq_along(y)
  
  df = do.call("rbind", lapply(time.labels, function(i){
    df           = data.frame(y[[i]]$AreaUncertainty)
    total_area   = sum(unlist(df$Mapped))
    df_m         = data.frame(t(df$Mapped))
    names(df_m)  = category.name
    df_m         = melt(df_m, measure.vars = names(df_m))
    df_m$Legend  = "Mapped"
    df_m$ci      = as.numeric(NA)
    df_m$Period  = i
    df_a         = data.frame(t(df$Adjusted))
    names(df_a)  = category.name
    df_a         = melt(df_a, measure.vars = names(df_a))
    df_a$Legend  = "Adjusted"
    df_a$ci      = as.numeric(df$ci)
    df_a$Period  = i
    df           = rbind(df_m, df_a)
    if(perc){
      df$value = df$value/total_area
      df$ci = df$ci/total_area
    }
    df
  }))
  
  limits = aes_string(ymax = "value + ci", ymin = "value - ci")
  dodge = position_dodge(width=0.9)
  
  gp = ggplot(df, aes_string(fill="Legend", y="value", x="variable")) + 
    facet_wrap(~Period, scales = "free") + 
    geom_bar(position="dodge", stat="identity", na.rm=TRUE) +
    geom_errorbar(limits, position=dodge, width=0.25, na.rm=TRUE) + 
    # scale_fill_grey(start = .6, end = .3) + 
    xlab("") + 
    ylab("Area")
  
  if(perc)
    gp = gp + scale_y_continuous(labels = percent)
  
  gp
  
}


