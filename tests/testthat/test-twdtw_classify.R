library(stars)
library(stringr)

samples <- st_read("inst/mato_grosso_brazil/samples.gpkg")

#tif_files <- system.file("inst/mato_grosso_brazil", package = "dtwSat") |>
#  dir(pattern = "\\.tif$", full.names = TRUE)
tif_files <- dir("inst/mato_grosso_brazil", pattern = "\\.tif$", full.names = TRUE)

# The acquisition date is in the file name are not the true acquisition date of each pixel
# MOD13Q1 is a 16-day composite product, so the acquisition date is the first day of the 16-day period
acquisition_date <- as.Date(str_extract(tif_files, "[0-9]{8}"), format = "%Y%m%d")

# Read the data as a stars object, setting time as a dimension and band as attribute
dc <- read_stars(tif_files, proxy = FALSE, along = list(time = acquisition_date)) |>
  st_set_dimensions(3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR", "DOY")) |>
  split("band")

# Remove the DOY band - this will be supported int the future
dc <- dc[c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")]

# Extract the samples
ts_samples <- st_extract(dc, samples)
ts_samples$id <- 1:dim(ts_samples)["geom"]
ts_samples$label <- samples$label

a <- as.data.frame(ts_samples)
a$geom <- NULL

head(a)

a[a$id == "1", ]

# training set ready!