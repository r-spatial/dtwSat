library(stars)
library(stringr)

# Read training samples
samples <- st_read(system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat"))

# Satellite image time sereis files
tif_files <- system.file("mato_grosso_brazil", package = "dtwSat") |>
  dir(pattern = "\\.tif$", full.names = TRUE)

# The acquisition date is in the file name are not the true acquisition date of each pixel
# MOD13Q1 is a 16-day composite product, so the acquisition date is the first day of the 16-day period
acquisition_date <- as.Date(str_extract(tif_files, "[0-9]{8}"), format = "%Y%m%d")

# Read the data as a stars object, setting time as a dimension and band as attribute
dc <- read_stars(tif_files, proxy = FALSE, along = list(time = acquisition_date)) |>
  st_set_dimensions(3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR", "DOY")) |>
  split("band")

# Remove the DOY band - this will be supported int the future
dc <- dc[c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")]

# Get temporal patters
ts_patterns <- create_patterns(x = dc, y = samples)

# Visualize patterns
plot_patterns(ts_patterns)
