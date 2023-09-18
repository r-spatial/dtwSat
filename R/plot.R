#' Plot Patterns from twdtw-knn1 model
#'
#' This function visualizes time series patterns from the \code{"twdtw_knn1"} model.
#' It produces a multi-faceted plot, where each facet represents a different time series
#' label from the model's data. Within each facet, different bands or indices (attributes)
#' are plotted as distinct lines, differentiated by color.
#'
#' @param x A model of class \code{"twdtw_knn1"}.
#'
#' @param bands A character vector specifying the bands or indices to plot.
#' If NULL (default), all available bands or indices in the data will be plotted.
#'
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggplot}}. Currently not used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object displaying the time series patterns.
#'
#' @seealso twdtw_knn1
#'
#' @inherit twdtw_knn1 examples
#'
#' @export
plot.twdtw_knn1 <- function(x, bands = NULL, ...) {

  # Convert the list of time series data into a long-format data.frame
  df <- x$data
  df$id <- 1:nrow(df)
  df <- unnest(df, cols = 'observations')

  # Select bands
  if(!is.null(bands)){
    df <- df[c('id', 'time', 'label', bands)]
  }

  # Pivote data into long format for ggplot2
  df <- pivot_longer(df, !c('id', 'label', 'time'), names_to = "band", values_to = "value")

  # Construct the ggplot
  gp <- ggplot(df, aes(x = .data$time, y = .data$value, colour = .data$band, group = interaction(.data$id, .data$band))) +
    geom_line() +
    facet_wrap(~label) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(title = "Bands")) +
    ylab("Value") +
    xlab("Time")

  return(gp)

}
