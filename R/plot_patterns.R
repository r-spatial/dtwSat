#' Plot Patterns from Time Series Data
#'
#' This function visualizes time series data patterns from the \code{"knn1_twdtw"} model.
#' It produces a multi-faceted plot,
#' where each facet represents a different time series pattern from the model's data.
#' Within each facet, different attributes (e.g., bands or indices) are plotted as
#' distinct lines, differentiated by color.
#'
#' @param x A model of class \code{"knn1_twdtw"}.
#' 
#' @param n An integer specifying the number of patterns to plot. Default is 12.
#' 
#' @param ... Additional arguments. Currently not used.
#' 
#' @details
#' The function is designed for visual inspection of time series patterns from the model. 
#' Each data frame in the model's data represents a time series. Columns within the data frame 
#' correspond to different attributes (e.g., bands or indices), while rows represent time points. 
#' The facet title is derived from the name of each time series in the model's data.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object displaying the time series patterns.
#'
#' @seealso knn1_twdtw
#'
#' @inherit knn1_twdtw examples
#'
#' @export
plot.knn1_twdtw <- function(x, n = 12, ...) {

  # Convert the list of time series data into a long-format data.frame
  df <- unnest(x$data[1:n, ], cols = 'observations')

  # Melt the data into long format suitable for ggplot2
  df <- pivot_longer(df, !c('label', 'time'), names_to = "band", values_to = "value")

  # Construct the ggplot
  gp <- ggplot(df, aes(x = .data$time, y = .data$value, colour = .data$band)) +
    geom_line() +
    facet_wrap(~label) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(title = "Bands")) +
    ylab("Value") +
    xlab("Time")

  return(gp)

}
