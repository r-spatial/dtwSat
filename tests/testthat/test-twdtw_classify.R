# Read training samples
samples <- st_read(system.file("mato_grosso_brazil/samples.gpkg", package = "dtwSat"), quiet = TRUE)

# Satellite image time sereis files
tif_files <- system.file("mato_grosso_brazil", package = "dtwSat") |>
  dir(pattern = "\\.tif$", full.names = TRUE)

# The acquisition date is in the file name are not the true acquisition date of each pixel
# MOD13Q1 is a 16-day composite product, so the acquisition date is the first day of the 16-day period
acquisition_date <- as.Date(regmatches(tif_files, regexpr("[0-9]{8}", tif_files)), format = "%Y%m%d")

# Read the data as a stars object setting the time/date for each observation using along
# This will prodcue a 4D array
# Split 'band' and 'time' dimensions to create a 2D array
dc <- read_stars(tif_files, proxy = FALSE, along = list(time = acquisition_date), RasterIO = list(bands = 1:6)) |>
  st_set_dimensions(3, c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")) |>
  split(c("band")) |>
  split(c("time"))

# Get temporal patters
m <- knn1_twdtw(x = dc, y = samples, formula = band ~ s(time), sampling_freq = 16)

# Visualize patterns
plot(m)

# Classify stars
system.time(lu <- predict(dc, model = m, drop_dimensions = TRUE, cycle_length = 'year', time_scale = 'day', time_weight = c(steepness = 0.1, midpoint = 50)))

# write_stars(lu, "~/Downloads/lu.tif")