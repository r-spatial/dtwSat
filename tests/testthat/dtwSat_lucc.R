###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2015-11-26                             #
#                                                             #
###############################################################


###############################################################
#### CREATE TEMPORAL PATTERNS AND APPLY TWDTW TO A STUDY AREA 


library(dtwSat)
library(raster)
library(rgdal)
library(parallel)
library(ggplot2)
library(caret)
library(lubridate)
library(scales)
library(reshape2)
library(Hmisc)


raster_dir = system.file("tif_brick_MT",  package="dtwSat")
ndvi  = brick(paste0(raster_dir, "/ndvi.tif"))
evi   = brick(paste0(raster_dir, "/evi.tif" ))
red   = brick(paste0(raster_dir, "/red.tif" ))
nir   = brick(paste0(raster_dir, "/nir.tif" ))
mir   = brick(paste0(raster_dir, "/mir.tif" ))
blue  = brick(paste0(raster_dir, "/blue.tif"))
doy   = brick(paste0(raster_dir, "/doy.tif"))
timeline = read.table(system.file("timeline.csv",  package="dtwSat"), as.is = TRUE)
timeline = as.Date(timeline[,1])
raster.list = createRasterTimeSeries(ndvi=ndvi, evi=evi, red=red, nir=nir, mir=mir, blue=blue, 
                             timeline=timeline, doy=doy)

### GET FILED SAMPLES 
field_samples = readOGR(system.file("samples_MT/samples_MT.shp",  package="dtwSat"), 
                        layer = "samples_MT")

table(field_samples$group)

### CREATE TEMPORAL PATTERNS 
# Extract time series from raster time series 
extract_time = system.time(
  ts.list <- extractTimeSeries(x = raster.list, y = field_samples, mc.cores = 4)
)

# I = row.names(field_samples)[field_samples$group=="Soybean-millet"]
# 
# df = do.call("rbind", lapply(I, function(i){
#   data.frame(time=index(ts.list[[i]]), ts.list[[i]])
# }))
# 
# ggplot(df, aes(time, evi)) + 
#   geom_point()


####
set.seed(1)
class_table.list = lapply(1:100, function(k){
  print(paste("Running simulation",k))
  # Split training and validation samples 
  I = createDataPartition(y = field_samples$group, p = 0.1, list = TRUE)$Resample1
  ts.training_sample   = ts.list[ I]
  ts.validation_sample = ts.list[-I]
  sp.training_sample   = field_samples[ I,]
  sp.validation_sample = field_samples[-I,]
  
  # Group land-use classes and shift time seres to the sam interval 
  groups = as.character(unique(sp.training_sample$group))
  names(groups) = groups
  J = lapply(groups, function(x) row.names(sp.training_sample)[sp.training_sample$group==x])
  
  # Create temporal patterns 
  samples.list = lapply(J, function(j) ts.training_sample[j] )
  patterns.list = lapply(samples.list, createPattern, from="2004-09-01", to="2005-09-01", 
                         freq=1, formula = y ~ s(time, bs="cc"))
  # grid.arrange(grobs=lapply(patterns.list, autoplot, facets = NULL), ncol=3)

  # Apply twdtw for the validation samples 
  validation_results = mclapply(ts.validation_sample, FUN=twdtw, patterns = patterns.list,
                                span = 180, weight.fun = logisticWeight(-0.1, 50), 
                                normalize.patterns = TRUE, patterns.length = 23,
                                mc.cores = 3)

  # Classify twdtw results
  s = names(validation_results)
  names(s) = s
  res = do.call("rbind", mclapply(s, function(i){
    pred = classifyIntervals(validation_results[[i]], from = sp.validation_sample[i,]$from,
                             to = sp.validation_sample[i,]$to, by = "12 month", 
                             overlap = 0.5, pattern.only = TRUE)
    data.frame(Reference = as.character(sp.validation_sample[i,]$group),
               Predicted = as.character(pred), stringsAsFactors = FALSE)
  }))
  return(res)
})


# Compute accuracy assessment 
assess_results = do.call("rbind", lapply(class_table.list, function(class_table){
  assess_table = table(Predicted=class_table$Predicted, Reference=class_table$Reference)
  user_accuracy = diag(assess_table) / rowSums(assess_table)
  prod_accuracy = diag(assess_table) / colSums(assess_table)
  data.frame(Group=names(user_accuracy), 
             User = user_accuracy,
             Producer = prod_accuracy)
}))

df = melt(assess_results, id="Group")
gp = ggplot(df, aes(x=Group, y=value)) + 
  # stat_summary(fun.data="median_hilow", width=0.5, geom="errorbar") + 
  stat_summary(fun.data="median_hilow", width=0.5, geom="crossbar", fill="grey") + 
  # stat_summary(fun.data="mean_sdl", width=0.5, geom="errorbar") + 
  # stat_summary(fun.data="mean_sdl", width=0.5, geom="crossbar", fill="grey") + 
  geom_point() +  
  # geom_boxplot() +
  coord_flip() + 
  xlab("") + 
  ylab("Accuracy %") + 
  facet_grid(. ~ variable) 
gp






# Run for the study area 
I = createDataPartition(y = field_samples$group, p = 0.1, list = TRUE)$Resample1
ts.training_sample   = ts.list[ I]
ts.validation_sample = ts.list[-I]
sp.training_sample   = field_samples[ I,]
sp.validation_sample = field_samples[-I,]

# Group land-use classes and shift time seres to the sam einterval 
groups = as.character(unique(sp.training_sample$group))
names(groups) = groups
J = lapply(groups, function(x) row.names(sp.training_sample)[sp.training_sample$group==x])

# Create temporal patterns 
samples.list = lapply(J, function(j) ts.training_sample[j] )
patterns.list = lapply(samples.list, createPattern, from="2004-09-01", to="2005-09-01", 
                       freq=1, formula = y ~ s(time, bs="cc"))

# Classification parameters 
twdtw_time = system.time(
    res <- twdtwApply( x=raster.list, patterns=patterns.list, 
                       from = "2007-09-01", to = "2013-09-01", by = "12 month", 
                       win.fun = classifyIntervals, win.size = c(1,1),
                       mc.cores = 3, chunk.size = 100, weight.fun = logisticWeight(-0.1,50), 
                       overlap = 0.5, normalize.patterns=TRUE, patterns.length = 23, 
                       pattern.only = TRUE)
)




layers = seq(from = as.Date("2007-09-01"), to = as.Date("2013-09-01"), by = "12 month")
layers = paste(format(layers[-length(layers)],"%b/%Y"), format(layers[-1]-1,"%b/%Y"), sep = "-")
res.df = data.frame(coordinates(res), res[] )
names(res.df) = c("Longitude", "Latitude", layers)
res.df = melt(res.df, id.vars = c("Longitude", "Latitude"))
res.df$value = factor(res.df$value, levels = c(seq_along(patterns.list), 255), labels = c(names(patterns.list), "Unclassified"))

class_colors = c("#996400", "#005500", "#d8b777", "#e6d219", "#e6bec8", "#c8c8c8")
names(class_colors) = c(names(patterns.list), "Unclassified")

gp = ggplot(data=res.df, aes_string(x="Longitude", y="Latitude")) +
  geom_raster(aes_string(fill="value")) + 
  scale_fill_manual(name="Land use", 
                    values = class_colors, guide = guide_legend(ncol=4)) + 
  facet_wrap(~variable, ncol=2) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(text = element_text(size = 8, family="Helvetica"),
        plot.margin=unit(c(1,1,0,0),"mm"),
        legend.position = "bottom") + 
  coord_fixed(ratio = 1)

gp


