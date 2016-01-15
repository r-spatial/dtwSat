## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  error = FALSE,
  results = "hide"
)

## ----study-area, echo = FALSE, eval = TRUE, fig.path='figure/', fig.width=6, fig.align='center'----
library(png)
library(grid)
grid.raster(readPNG("study_area.png"))

## ---- echo = TRUE, eval = TRUE-------------------------------------------
library(dtwSat)
library(raster)
raster_dir = system.file('lucc_MT',  package = 'dtwSat')
raster_files_list = paste(raster_dir, dir(raster_dir, pattern = '.tif$'), sep = '/')

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  bands = 1:6
#  names(bands) = c('ndvi','evi','red','nir','mir','blue')
#  stack_list = lapply(bands, function(x)
#            stack(raster_files_list, bands = x)
#    )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  doy = stack(raster_files_list, bands = 7)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
timeline = read.table(system.file('lucc_MT/timeline.csv', package = 'dtwSat'), as.is = TRUE)

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  raster_timeseries = createRasterTimeSeries(x = stack_list, timeline = timeline[,1], doy = doy)

## ---- echo = FALSE, eval = TRUE------------------------------------------
# Fast version using raster brick  
flist = dir("~/lucc_brick/", pattern = ".tif$")
fnames = unlist(lapply(flist, function(x) unlist(strsplit(x, split="[.]"))[1]))
flist = paste0("~/lucc_brick/",flist)
names(flist) = fnames
brick_list = lapply(flist[-2], brick)
doy = brick(flist[2])
raster_timeseries = createRasterTimeSeries(x = brick_list, timeline = timeline[,1], doy = doy)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
library(rgdal)
field_samples = readOGR(system.file('lucc_MT/samples_MT.shp', package="dtwSat"),
                        layer = "samples_MT")

## ---- echo = TRUE, eval = TRUE-------------------------------------------
names(field_samples)
table(field_samples$group)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
ts_list = extractTimeSeries(x = raster_timeseries, y = field_samples, mc.cores = 1)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
groups = as.character(unique(field_samples$group))
names(groups) = groups
field_samples_list = lapply(groups, function(x) ts_list[field_samples$group==x])
patterns_list = lapply(field_samples_list, createPattern, 
                       freq = 8, formula = y ~ s(time, bs = "cc"))

## ----temporal-patterns, echo = TRUE, eval = TRUE, fig.path='figure/', fig.width=6, fig.height=4, fig.align='center'----
library(ggplot2)
plotPatterns(patterns_list) + 
  theme(text = element_text(size = 8, family="Helvetica"),
        legend.position = "bottom")

## ---- echo = TRUE, eval = TRUE-------------------------------------------
## Time-weight function for TWDTW analysis. See ?twdtw for details 
weight.fun = logisticWeight(alpha=-0.1, beta=50)
## Classification function. See ?classifyIntervals for details
win.fun = classifyIntervals 
## Classification intervals and overlap. See ?classifyIntervals for details
breaks = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), by = "12 month")
overlap = 0.5 
# Legend parameter. See ?classifyIntervals for details
levels = c(seq_along(patterns_list), 255)
labels = c(names(patterns_list), "Unclassified")
colors = c("#996400", "#005500", "#D8B777", "#E6D219", "#E6BEC8", "#C8C8C8")
names(colors) = labels

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  land_use_maps = twdtwApply(x = raster_timeseries, patterns = patterns_list,
#                             mc.cores = 4, win.fun = win.fun, weight.fun = weight.fun,
#                             breaks = breaks, overlap = overlap, levels = levels,
#                             labels = labels, simplify= TRUE)

## ---- echo = FALSE, eval = FALSE-----------------------------------------
#  proctime = system.time(
#     land_use_maps <- twdtwApply(x = raster_timeseries, patterns = patterns_list,
#                                 mc.cores = 4, win.fun = win.fun, weight.fun = weight.fun,
#                                 breaks = breaks, overlap = overlap, levels = levels,
#                                 labels = labels, simplify= TRUE)
#  )
#  proctime/60
#  #writeRaster(land_use_maps, "./inst/lucc_MT/lucc_MT.tif", overwrite=TRUE)

## ---- echo = FALSE, eval = TRUE------------------------------------------
### Load classification not to run again 
land_use_maps = brick(system.file('lucc_MT/lucc_MT.tif', package = 'dtwSat'))

## ----plot-map, echo = TRUE, eval = TRUE, fig.path='figure/', fig.width=6, fig.height=4, fig.align='center'----
plotLUCC(x = land_use_maps, type = "map", layer.labels = 2008:2013, 
         levels = levels, labels = labels, colors = colors) + 
         theme(text = element_text(size = 8, family="Helvetica"))

## ----plot-area, echo = TRUE, eval = TRUE, fig.path='figure/', fig.width=6, fig.height=4, fig.align='center'----
plotLUCC(x = land_use_maps, type = "area", layer.labels = 2008:2013, 
         levels = levels, labels = labels, 
         colors = colors) + 
         theme(text = element_text(size = 8, family="Helvetica"))

## ----plot-change, echo = TRUE, eval = TRUE, fig.path='figure/', fig.width=6, fig.height=4, fig.align='center'----
plotLUCC(x = land_use_maps, type = "change", layer.labels = 2008:2013, 
         levels = levels, labels = labels, colors = colors) + 
         theme(text = element_text(size = 8, family="Helvetica"))

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  library(caret)
#  library(parallel)
#  library(reshape2)
#  n = 100 # Number of repetitions
#  p = 0.1 # p% Training (1-p)% Validation
#  set.seed(1)
#  training_sample = createDataPartition(y = field_samples$group,
#                                        times = n, p = 0.1, list = TRUE)
#  
#  assess_list = lapply(training_sample, function(I){
#    cat(".")
#    # Split training and validation samples
#    ts_training_sample   = ts_list[ I]
#    ts_validation_sample = ts_list[-I]
#    sp_training_sample   = field_samples[ I,]
#    sp_validation_sample = field_samples[-I,]
#  
#    # Group land-use classes
#    groups = as.character(unique(sp_training_sample$group))
#    names(groups) = groups
#    J = lapply(groups, function(x)
#      row.names(sp_training_sample)[sp_training_sample$group==x])
#  
#    # Create temporal patterns
#    samples_list = lapply(J, function(j) ts_training_sample[j] )
#    patterns_list = lapply(samples_list, createPattern, freq=8,
#                           from="2007-09-01", to="2008-09-01",
#                           formula = y ~ s(time, bs="cc"))
#    # grid.arrange(grobs=lapply(patterns.list, autoplot, facets = NULL), ncol=3)
#  
#    # Apply twdtw for the validation samples
#    validation_results = mclapply(ts_validation_sample, FUN=twdtw,
#                                patterns = patterns_list,
#                                weight.fun = weight.fun, mc.cores=3)
#  
#    # Classify twdtw results
#    s = names(validation_results)
#    names(s) = s
#    res = do.call("rbind", mclapply(s, mc.cores=3, function(i){
#      from = sp_validation_sample[i,]$from
#      to = sp_validation_sample[i,]$to
#      pred = classifyIntervals(validation_results[[i]],
#                               from = from, to = to,
#                               by = "12 month",
#                               overlap = overlap)
#      data.frame(Reference = as.character(sp_validation_sample[i,]$group),
#                 Predicted = as.character(pred$pattern), stringsAsFactors = FALSE)
#    }))
#    return(res)
#  })

## ---- echo = FALSE, eval = FALSE-----------------------------------------
#  save(assess_list, file="./inst/lucc_MT/assess_list.RData")

## ---- echo = FALSE, eval = TRUE------------------------------------------
## Load classification not to run again 
load(system.file('lucc_MT/assess_list.RData', package = 'dtwSat'))

## ---- echo = TRUE, eval = TRUE-------------------------------------------
assess_results = do.call("rbind", lapply(assess_list, function(class_table){
  assess_table = table(Predicted=class_table$Predicted, Reference=class_table$Reference)
  user_accuracy = diag(assess_table) / rowSums(assess_table)
  prod_accuracy = diag(assess_table) / colSums(assess_table)
  data.frame(Group=names(user_accuracy), 
             User = user_accuracy,
             Producer = prod_accuracy)
}))

## ----plot-accuracy, echo = TRUE, eval = TRUE, fig.path='figure/', fig.width=6, fig.height=4, fig.align='center'----
library(reshape2)
library(scales)
df = melt(assess_results, id="Group")
ggplot(df, aes(x=Group, y=value)) + 
  stat_summary(fun.data="median_hilow", width=0.5, geom="crossbar", fill="grey") + 
  geom_point() +  
  facet_grid(. ~ variable) + 
  scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
  xlab("") + 
  ylab("Accuracy") + 
  coord_flip()

