library(stars)
library(stringr)
library(dplyr)

sf_use_s2(FALSE)

# Data extracted from MOD13Q1 tile h12v10

sd <- sf::gdal_subdatasets(hdf_file)


samples <- st_read("inst/mato_grosso_brazil/samples.gpkg")
table(samples$label)

modis_dates <- scan("inst/lucc_MT/data/timeline", what = "character")

evi  <- raster::brick(system.file("lucc_MT/data/evi.tif",  package = "dtwSat"))
ndvi <- raster::brick(system.file("lucc_MT/data/ndvi.tif", package = "dtwSat"))
red  <- raster::brick(system.file("lucc_MT/data/red.tif",  package = "dtwSat"))
blue <- raster::brick(system.file("lucc_MT/data/blue.tif", package = "dtwSat"))
nir  <- raster::brick(system.file("lucc_MT/data/nir.tif",  package = "dtwSat"))
mir  <- raster::brick(system.file("lucc_MT/data/mir.tif",  package = "dtwSat"))
doy  <- raster::brick(system.file("lucc_MT/data/doy.tif",  package = "dtwSat"))

# filter modis dates between 2011-09-01 and 2012-08-31
for (i in which(modis_dates >= "2011-09-01" & modis_dates <= "2012-08-31")){
  r <- stack(subset(evi, i), subset(ndvi, i), subset(red, i), subset(blue, i), subset(nir, i), subset(mir, i), subset(doy, i))
  names(r) <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR", "DOY")
  writeRaster(r, filename = paste0("inst/mato_grosso_brazil/MOD13Q1_", str_remove_all(modis_dates[i], "-"), "_subset_from_h12v10.tif"),
              overwrite = TRUE, wopt= list(gdal=c("COMPRESS=DEFLATE", "TFW=YES")))
}
