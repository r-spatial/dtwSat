library(mgcv)

# Read training samples
samples <- st_read(system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat"), quiet = TRUE)

# Satellite image time sereis files
tif_files <- dir(system.file("mato_grosso_brazil", package = "dtwSat"), pattern = "\\.tif$", full.names = TRUE)

# The acquisition date is in the file name are not the true acquisition date
# of each pixel. MOD13Q1 is a 16-day composite product, so the acquisition date
# is the first day of the 16-day period
acquisition_date <- as.Date(regmatches(tif_files, regexpr("[0-9]{8}", tif_files)), format = "%Y%m%d")

# Read the data as a stars object setting the time/date for each observation
# using along. This will prodcue a 4D array (data-cube) which will then be converted
# to a 3D array by spliting the 'band' dimension
dc <- read_stars(tif_files,
                 proxy = FALSE,
                 along = list(time = acquisition_date),
                 RasterIO = list(bands = 1:6))

dc <- st_set_dimensions(dc, 3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR"))
dc <- split(dc, c("band"))

# Create a knn1-twdtw model
system.time(
  m <- twdtw_knn1(x = dc,
                  y = samples,
                  smooth_fun = function(x, y) gam(y ~ s(x), data = data.frame(x = x, y = y)),
                  cycle_length = 'year',
                  time_scale = 'day',
                  time_weight = c(steepness = 0.1, midpoint = 50))
)

print(m)

# Visualize model patterns
plot(m)

# Classify satellite image time series
system.time(lu <- predict(dc, model = m))

print(lu)

# # Visualise land use classification
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

# ### OTHER TESTS
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = function(x, y) gam(y ~ s(x), data = data.frame(x = x, y = y)),
                resampling_freq = 60)

plot(m)

system.time(
  lu <- predict(dc,
                model = m,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50))
)

print(lu)

# Visualise land use classification
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

# Test model without samples reduction
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = NULL)

print(m)

plot(m)

plot(m, bands = c('EVI', 'NDVI'))

# Test custom smooth function -- mean value for each time in dc
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = function(x, y) lm(y ~ factor(x), data = data.frame(x = x, y = y)))

print(m)

plot(m)

plot(m, bands = c('EVI', 'NDVI'))
