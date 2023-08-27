#' Plot Patterns from Time Series Data
#'
#' This function takes a list of time series data and creates a multi-faceted plot
#' where each facet corresponds to a different time series from the list. 
#' Within each facet, different attributes (columns of the time series) are 
#' plotted as lines with different colors.
#'
#' @param x A list where each element is a data.frame representing a time series.
#'          Each data.frame should have the same number of rows and columns, 
#'          with columns representing different attributes (e.g., bands or indices)
#'          and rows representing time points.
#'          The name of each element in the list will be used as the facet title.
#' 
#' @param ... Not used.
#' 
#' @return A ggplot object displaying the time series patterns.
#'
#' @export
plot_patterns = function(x, ...) {
   
  # Convert the list of time series data into a long-format data.frame
  df.p = do.call("rbind", lapply(names(x), function(p) {
    ts = x[[p]]
    # Create a new data.frame with a 'Time' column and a 'Pattern' column
    # representing the name of the current time series (facet name).
    data.frame(Time = 1:nrow(ts), ts, Pattern = p) # Assuming the time series are evenly spaced
  }))
  
  # Melt the data into long format suitable for ggplot2
  df.p = melt(df.p, id.vars = c("Time", "Pattern"))
  
  # Construct the ggplot
  gp = ggplot(df.p, aes(x = .data$Time, y = .data$value, colour = .data$variable)) + 
    geom_line() + 
    facet_wrap(~Pattern) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(title = "Bands")) + 
    ylab("Value")
  
  return(gp)
}