cat("Loading mgcv...\n")
library(mgcv)

cat("Reading training samples...\n")
samples <- st_read(system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat"))

cat("Geting satellite image time series files...\n")
tif_files <- dir(system.file("mato_grosso_brazil", package = "dtwSat"), pattern = "\\.tif$", full.names = TRUE)

cat("Geting acquisition dates from file names...\n")
acquisition_date <- as.Date(regmatches(tif_files, regexpr("[0-9]{8}", tif_files)), format = "%Y%m%d")

cat("Reading data as stars object...\n")
dc <- read_stars(tif_files,
                 proxy = FALSE,
                 along = list(time = acquisition_date),
                 RasterIO = list(bands = 1:6))

cat("Setting band names...\n")
dc <- st_set_dimensions(dc, 3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR"))

cat("Converting 4D into a 3D datacube...\n")
dc <- split(dc, c("band"))

cat("Creating twdtw_knn1 model using GAM...\n")
system.time(
  m <- twdtw_knn1(x = dc,
                  y = samples,
                  smooth_fun = function(x, y) gam(y ~ s(x), data = data.frame(x = x, y = y)),
                  cycle_length = 'year',
                  time_scale = 'day',
                  time_weight = c(steepness = 0.1, midpoint = 50))
)

cat("Printing model...\n")
print(m)

cat("Visualizing model patterns...\n")
plot(m)

cat("Classifying satellite image time series...\n")
system.time(lu <- predict(dc, model = m))

cat("Printing land use classification...\n")
print(lu)

cat("Visualizing land use classification...\n")
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

cat("Testing model with resampling frequency of 60 and GAM...\n")
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = function(x, y) gam(y ~ s(x), data = data.frame(x = x, y = y)),
                resampling_freq = 60)

cat("Visualizing resampled patterns...\n")
plot(m)

cat("Classifying satellite image time series with resampled patterns...\n")
system.time(
  lu <- predict(dc,
                model = m,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50))
)

cat("Printing land use classification of resampled patterns...\n")
print(lu)

cat("Visualizing land use classification of resampled patterns...\n")
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

cat("Testing model without sample reduction...\n")
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = NULL)

cat("Printing model without sample reduction...\n")
print(m)

cat("Visualizing model patterns without sample reduction...\n")
plot(m)

cat("Visualizing model patterns without sample reduction for EVI and NDVI...\n")
plot(m, bands = c('EVI', 'NDVI'))

cat("Testing custom smooth function using the average for each observation date...\n")
m <- twdtw_knn1(x = dc,
                y = samples,
                cycle_length = 'year',
                time_scale = 'day',
                time_weight = c(steepness = 0.1, midpoint = 50),
                smooth_fun = function(x, y) lm(y ~ factor(x), data = data.frame(x = x, y = y)))

cat("Printing model with custom smooth...\n")
print(m)

cat("Visualizing model patterns with custom smooth...\n")
plot(m)

cat("Visualizing model patterns with custom smooth for EVI and NDVI...\n")
plot(m, bands = c('EVI', 'NDVI'))

cat("Classifying satellite image time series with custom smooth...\n")
system.time(lu <- predict(dc, model = m))

cat("Printing land use classification with custom smooth...\n")
print(lu)

cat("Visualizing land use classification with custom smooth...\n")
ggplot() +
  geom_stars(data = lu) +
  theme_minimal()

cat("Script completed.\n")
