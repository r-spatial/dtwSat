% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.R
\name{twdtw_knn1}
\alias{twdtw_knn1}
\title{Train a KNN-1 TWDTW model with optional GAM resampling}
\usage{
twdtw_knn1(
  x,
  y,
  time_weight,
  cycle_length,
  time_scale,
  formula = NULL,
  start_column = "start_date",
  end_column = "end_date",
  label_colum = "label",
  sampling_freq = NULL,
  ...
)
}
\arguments{
\item{x}{A three-dimensional stars object (x, y, time) with bands as attributes.}

\item{y}{An sf object with the coordinates of the training points.}

\item{time_weight}{A numeric vector with length two (steepness and midpoint of logistic weight) or a function.
See details in \link[twdtw]{twdtw}.}

\item{cycle_length}{The length of the cycle, e.g. phenological cycles. Details in \link[twdtw]{twdtw}.}

\item{time_scale}{Specifies the time scale for the observations. Details in \link[twdtw]{twdtw}.}

\item{formula}{Either NULL or a formula to reduce samples of the same label using Generalized Additive Models (GAM).
Default is \code{band ~ s(time)}. See details.}

\item{start_column}{Name of the column in y that indicates the start date. Default is 'start_date'.}

\item{end_column}{Name of the column in y that indicates the end date. Default is 'end_date'.}

\item{label_colum}{Name of the column in y containing land use labels. Default is 'label'.}

\item{sampling_freq}{The time frequency for sampling, including the unit (e.g., '16 day').
If NULL, the function will infer the frequency. This parameter is only used if a formula is provided.}

\item{...}{Additional arguments passed to the \link[mgcv]{gam} function and to \link[twdtw]{twdtw} function.}
}
\value{
A 'twdtw_knn1' model containing the trained model information and the data used.
}
\description{
This function prepares a KNN-1 model with the Time Warp Dynamic Time Warping (TWDTW) algorithm.
If a formula is provided, the training samples are resampled using Generalized Additive Models (GAM).
}
\details{
If \code{formula} is NULL, the KNN-1 model will retain all training samples. If a formula is passed (e.g., \code{band ~ \link[mgcv]{s}(time)}),
then samples of the same label (land cover class) will be resampled using GAM.
Resampling can significantly reduce prediction processing time.
}
\examples{
\dontrun{

# Read training samples
samples_path <- 
  system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat")

samples <- st_read(samples_path, quiet = TRUE)

# Get satellite image time sereis files
tif_path <- system.file("mato_grosso_brazil", package = "dtwSat")
tif_files <- dir(tif_path, pattern = "\\\\.tif$", full.names = TRUE)

# Get acquisition dates
acquisition_date <- regmatches(tif_files, regexpr("[0-9]{8}", tif_files))
acquisition_date <- as.Date(acquisition_date, format = "\%Y\%m\%d")

# Create a 3D datacube
dc <- read_stars(tif_files,
                 proxy = FALSE,
                 along = list(time = acquisition_date),
                 RasterIO = list(bands = 1:6))
dc <- st_set_dimensions(dc, 3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR"))
dc <- split(dc, c("band"))

# Create a knn1-twdtw model
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                formula = band ~ s(time))

print(m)

# Visualize model patterns
plot(m)

# Classify satellite images
system.time(lu <- predict(dc, model = m))

# Visualise land use classification
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

}
}
