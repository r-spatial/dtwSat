# ###############################################################
# #                                                             #
# #   (c) Victor Maus <vwmaus1@gmail.com>                       #
# #       Institute for Geoinformatics (IFGI)                   #
# #       University of Muenster (WWU), Germany                 #
# #                                                             #
# #       Earth System Science Center (CCST)                    #
# #       National Institute for Space Research (INPE), Brazil  #
# #                                                             #
# #                                                             #
# #   R Package dtwSat - 2015-11-26                             #
# #                                                             #
# ###############################################################
# 
# 
# ###############################################################
# #### CREATE TEMPORAL PATTERNS AND APPLY TWDTW TO A STUDY AREA 
# 
# 
# library(dtwSat)
# library(raster)
# library(rgdal)
# library(parallel)
# library(ggplot2)
# library(caret)
# library(lubridate)
# library(scales)
# library(reshape2)
# 
# 
# raster_dir = system.file("tif_MT",  package="dtwSat")
# ndvi  = brick(paste0(raster_dir, "/ndvi.tif"))
# evi   = brick(paste0(raster_dir, "/evi.tif" ))
# red   = brick(paste0(raster_dir, "/red.tif" ))
# nir   = brick(paste0(raster_dir, "/nir.tif" ))
# mir   = brick(paste0(raster_dir, "/mir.tif" ))
# blue  = brick(paste0(raster_dir, "/blue.tif"))
# doy   = brick(paste0(raster_dir, "/doy.tif"))
# timeline = read.table(system.file("tif_MT/timeline.csv",  package="dtwSat"), as.is = TRUE)
# timeline = as.Date(timeline[,1])
# raster.list = createRasterTS(ndvi=ndvi, evi=evi, red=red, nir=nir, mir=mir, blue=blue, 
#                              timeline=timeline, doy=doy)
# 
# ### GET FILED SAMPLES 
# field_samples = readOGR(system.file("samples_MT/samples_MT.shp",  package="dtwSat"), 
#                         layer = "samples_MT")
# 
# 
# ### CREATE TEMPORAL PATTERNS 
# # Extract time series from raster time series 
# extract_time = system.time(
#   ts.list <- extractSampleTimeSeries(x = raster.list, y = field_samples, mc.cores = 4)
# )
# 
# ####  ####### MAKE A LOOP MANY ITERACTIONS KEEPING ACCURACY RESULTS 
# 
# ### SPLIT TRAINING AND VALIDATION SAMPLES 
# 
# set.seed(1)
# I = createDataPartition(y = field_samples$group, p = 0.1, list = TRUE)$Resample1
# table(field_samples$group[I])
# ts.training_sample   = ts.list[ I]
# ts.validation_sample = ts.list[-I]
# sp.training_sample   = field_samples[ I,]
# sp.validation_sample = field_samples[-I,]
# 
# # Group land-use classes and shift time seres to the sam einterval 
# groups = as.character(unique(sp.training_sample$group))
# names(groups) = groups
# J = lapply(groups, function(x) row.names(sp.training_sample)[sp.training_sample$group==x])
# 
# # Create temporal patterns 
# s = lapply(J, function(j) ts.training_sample[j] )
# patterns.list = lapply(s, createPattern, from="2004-09-01", to="2005-09-01", freq=1, 
#                        formula = y ~ s(time, bs="cc"))
# grid.arrange(grobs=lapply(patterns.list, autoplot, facets = NULL), ncol=3)
# 
# 
# # Classification parameters 
# twdtw_time = system.time(
#     res <- twdtwApply( x=raster.list, patterns=patterns.list, 
#                        from = "2007-09-01", to = "2013-09-01", by = "12 month", 
#                        win.fun = classifyIntervals, win.size = c(1,1),
#                        mc.cores = 3, chunk.size = 100, weight.fun = logisticWeight(-0.1,100), 
#                        overlap = 0.5, normalize.patterns=TRUE, patterns.length = 23, 
#                        pattern.only = TRUE)
# )
# 
# 
# 
# 
# # Map colours 
# group_colors = rgb(class_table$R, class_table$G, class_table$B, maxColorValue = 255)
# 
# layers = seq(from = as.Date("2007-09-01"), to = as.Date("2013-09-01"), by = "12 month")
# layers = paste(format(layers[-length(layers)],"%b/%Y"), format(layers[-1]-1,"%b/%Y"), sep = "-")
# res.df = data.frame(coordinates(res), res[] )
# names(res.df) = c("Longitude", "Latitude", layers)
# res.df = melt(res.df, id.vars = c("Longitude", "Latitude"))
# res.df$value = factor(res.df$value, levels = c(seq_along(patterns.list), 255), labels = c(names(patterns.list), "Unclassified"))
# 
# class_colors = c("#996400", "#005500", "#d8b777", "#e6d219", "#e6bec8", "#c8c8c8")
# names(class_colors) = c(names(patterns.list), "Unclassified")
# 
# gp = ggplot(data=res.df, aes_string(x="Longitude", y="Latitude")) +
#   geom_raster(aes_string(fill="value")) + 
#   scale_fill_manual(name="Land use", 
#                     values = class_colors, guide = guide_legend(ncol=4)) + 
#   facet_wrap(~variable, ncol=2) + 
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) + 
#   theme(text = element_text(size = 8, family="Helvetica"),
#         plot.margin=unit(c(1,1,0,0),"mm"),
#         legend.position = "bottom") + 
#   coord_fixed(ratio = 1)
# 
# gp
# 
