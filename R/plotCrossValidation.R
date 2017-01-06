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
#   R Package dtwSat - 2016-11-29                             #
#                                                             #
###############################################################


#' @title Plotting cross-validation 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description Method for plotting cross-validation results.
#' 
#' @param x An object of class \code{\link[dtwSat]{plotCrossValidation}}.
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwCrossValidation}}
#'  
#' @examples
#' \dontrun{
#' # Data folder 
#' data_folder = system.file("lucc_MT/data", package = "dtwSat")
#' 
#' # Read dates 
#' dates = scan(paste(data_folder,"timeline", sep = "/"), what = "dates")
#' 
#' # Read raster time series 
#' evi = brick(paste(data_folder,"evi.tif", sep = "/"))
#' raster_timeseries = twdtwRaster(evi, timeline = dates)
#' 
#' # Read field samples 
#' field_samples = read.csv(paste(data_folder,"samples.csv", sep = "/")) 
#' table(field_samples[["label"]])
#' 
#' # Read field samples projection 
#' proj_str = scan(paste(data_folder,"samples_projection", sep = "/"), 
#'      what = "character")
#' 
#' # Get sample time series from raster time series 
#' field_samples_ts = getTimeSeries(raster_timeseries, 
#'      y = field_samples, proj4string = proj_str)
#' field_samples_ts
#' 
#' # Run cross validation
#' set.seed(1)
#' # Define TWDTW weight function 
#' log_fun = logisticWeight(alpha=-0.1, beta=50) 
#' cross_validation = twdtwCrossValidation(field_samples_ts, times=3, p=0.1, 
#'                           freq = 8, formula = y ~ s(x, bs="cc"), weight.fun = log_fun)
#' 
#' summary(cross_validation)
#' 
#' plot(cross_validation)
#' 
#' }
#' 
#' @export
plotCrossValidation = function(x){
  UA = do.call("rbind", lapply(x@accuracy, function(x) data.frame(label="UA", rbind(x$UsersAccuracy))))
  PA = do.call("rbind", lapply(x@accuracy, function(x) data.frame(label="PA", rbind(x$UsersAccuracy))))
  df = melt(rbind(UA,PA), id="label")
  df$label = factor(df$label, levels = c("UA", "PA"), 
                       labels = c("User's Accuracy", "Producer's Accuracy"))
  df$variable = factor(df$variable, levels = levels(df$variable), 
                    labels = gsub("[.]","-",levels(df$variable)))
  gp = ggplot(df, aes_string(x="variable", y="value")) +
        stat_summary(fun.data="mean_cl_boot", fun.args=list(conf.int = .99), 
                     width=0.5, geom="crossbar", size=0.1, fill = "gray") + 
        geom_point(size=0.2) +  
        facet_grid(. ~ label) + 
        scale_y_continuous(limits = c(0,1), labels = percent, breaks = seq(0,1,.2)) + 
        xlab("") + 
        ylab("Accuracy") + 
        coord_flip()
  gp
}


