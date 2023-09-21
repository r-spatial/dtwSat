# Read training samples
samples <- system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat") |>
  st_read(quiet = TRUE)

# Satellite image time sereis files
tif_files <- system.file("mato_grosso_brazil", package = "dtwSat") |>
  dir(pattern = "\\.tif$", full.names = TRUE)

# The acquisition date is in the file name are not the true acquisition date
# of each pixel. MOD13Q1 is a 16-day composite product, so the acquisition date
# is the first day of the 16-day period
acquisition_date <- regmatches(tif_files, regexpr("[0-9]{8}", tif_files)) |>
  as.Date(format = "%Y%m%d")

# Read the data as a stars object setting the time/date for each observation
# using along. This will prodcue a 4D array (data-cube) which will then be converted
# to a 3D array by spliting the 'band' dimension
dc <- read_stars(tif_files,
                 proxy = FALSE,
                 along = list(time = acquisition_date),
                 RasterIO = list(bands = 1:6)) |>
  st_set_dimensions(3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")) |>
  split(c("band"))

# Create a knn1-twdtw model
system.time(
  m <- twdtw_knn1(x = dc,
                  y = samples,
                  cycle_length = 'year',
                  time_scale = 'day',
                  time_weight = c(steepness = 0.1, midpoint = 50),
                  formula = band ~ s(time))
)

print(m)

# Visualize model patterns
plot(m)

# Classify satellite image time series
system.time(lu <- predict(dc, model = m))

# Visualise land use classification
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

### OTHER TESTS
# split time first
dc <- read_stars(tif_files,
                 proxy = FALSE,
                 along = list(time = acquisition_date),
                 RasterIO = list(bands = 1:6)) |>
  st_set_dimensions(3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")) |>
  split(c("band"))

m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                formula = band ~ s(time),
                sampling_freq = 60)

plot(m)

system.time(
  lu <- predict(dc,
                model = m,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50))
)

# Visualise land use classification
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

# Test model without samples reduction
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50))

plot(m)

plot(m, bands = c('EVI', 'NDVI'))
