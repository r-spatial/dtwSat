# Example of TWDTW analysis using time series extracted from raster 
library(dtwSat)
library(caret) 

# Load raster data 
evi  <- brick(system.file("lucc_MT/data/evi.tif",  package = "dtwSat"))
ndvi <- brick(system.file("lucc_MT/data/ndvi.tif", package = "dtwSat"))
red  <- brick(system.file("lucc_MT/data/red.tif",  package = "dtwSat"))
blue <- brick(system.file("lucc_MT/data/blue.tif", package = "dtwSat"))
nir  <- brick(system.file("lucc_MT/data/nir.tif",  package = "dtwSat"))
mir  <- brick(system.file("lucc_MT/data/mir.tif",  package = "dtwSat"))
doy  <- brick(system.file("lucc_MT/data/doy.tif",  package = "dtwSat"))
timeline <- 
  scan(system.file("lucc_MT/data/timeline", package = "dtwSat"), what="date")

# Create raster time series 
rts <- twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)

# Load field samples and projection 
field_samples <- 
  read.csv(system.file("lucc_MT/data/samples.csv", package = "dtwSat"))
proj_str <- 
  scan(system.file("lucc_MT/data/samples_projection", package = "dtwSat"), 
       what = "character")

# Split samples for training (10%) and validation (90%) using stratified sampling 
set.seed(1)
I <- unlist(createDataPartition(field_samples$label, p = 0.1))
training_samples <- field_samples[I, ]
validation_samples <- field_samples[-I, ]

# Get time series form raster
training_ts <- getTimeSeries(rts, y = training_samples, proj4string = proj_str)
validation_ts <- getTimeSeries(rts, y = validation_samples, proj4string = proj_str)

# Create temporal patterns 
temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))

# Set TWDTW weight function 
log_fun <- logisticWeight(-0.1, 50)

# Run serial TWDTW analysis 
r_twdtw <- twdtwApply(x = validation_ts, 
                     y = temporal_patterns, weight.fun = log_fun)

# Classify raster based on the TWDTW analysis 
r_lucc <- twdtwClassify(x = r_twdtw, overlap = 0.5)

# Accuracy assessment 
twdtw_assess <- twdtwAssess(r_lucc, area = 53664.67, conf.int = .95)

# Plot map accuracy 
plot(twdtw_assess, type = "accuracy")

# Plot area uncertainty 
plot(twdtw_assess, type = "area")

# Get latex table with error matrix 
twdtwXtable(twdtw_assess, table.type = "matrix")

# Get latex table with error accuracy 
twdtwXtable(twdtw_assess, table.type = "accuracy")

# Get latex table with area uncertainty 
twdtwXtable(twdtw_assess, table.type = "area")

