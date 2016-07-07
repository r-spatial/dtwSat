## ---- echo=FALSE, cache = FALSE------------------------------------------
library(knitr)
opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  error = FALSE,
  results = "hide",
  cache = TRUE,
  comment = ""
)
#knit_hooks$set(purl = hook_purl) # Save chunks to Rscript file  

## ---- echo=FALSE, eval = TRUE, cache = FALSE-----------------------------
# Install dtwSat package  
#install.packages("dtwSat_0.2.0.9000.tar.gz")

## ---- echo=FALSE, eval = TRUE, cache = FALSE-----------------------------
library(dtwSat)
library(ggplot2)
library(scales)
library(reshape2)

new_theme = theme_get()
new_theme$text$family = "Helvetica"
new_theme$text$size = 8
old_theme = theme_set(new_theme)

tab_format = "latex"
page_width = 5.590551#in 14.2#cm
page_height = 9.173228#in 23.3#cm

## ----twdtw-example, echo = FALSE, eval = TRUE, fig.width=page_width, fig.height=page_height/3.5, fig.align='center', fig.cap='Matches of the known temporal pattern to subintervals of the long-term time series. The solid black line is the long-term time series, the colored lines are the different matches of the same pattern ordered by TWDTW dissimilarity measure, and the gray dashed lines are the matching points.', fig.pos='h'----
n=4
log_fun = logisticWeight(alpha = -0.1, beta = 100)
ts = twdtwTimeSeries(example_ts, labels="Time series")
patt = twdtwTimeSeries(patterns.list$Cotton, labels="Matches")
mat = twdtwApply(x=ts, y=patt, weight.fun=log_fun, keep=TRUE, n=n)
df_dist = mat[[1]]
df_dist$label = paste("Distance:",round(df_dist$distance,2))
df_dist$y = 1.8 
plotMatches(mat, attr="evi", k=n) + 
  ylab("Time series                  Pattern") + 
  geom_text(data=df_dist, mapping = aes_string(x='to', y='y', label='label'),
            size = 2, family="Helvetica") + 
  theme(legend.position="none")

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
ts = twdtwTimeSeries(example_ts, labels="Time series")
patterns_ts = twdtwTimeSeries(patterns.list)
example_ts_labels

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
library(dtwSat)
names(patterns.list)
head(example_ts, n = 2)

## ----example-timeseries, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_height/3, fig.align='center', fig.cap='Example of time series based on MODIS product MOD13Q1 \\citep{Friedl:2010}. The labels of the phenological cycle are shown in the plot.', fig.pos='!h'----
plot(ts, type = "timeseries") + 
  annotate(geom = "text", x = example_ts_labels$from+90, y = 0.98, 
  label = example_ts_labels$label, size = 2)

## ----temporal-patterns-soy-cot-mai, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_height/3.5, fig.align='center', fig.cap='Temporal patterns of soybean, cotton, and maize based on MODIS product MOD13Q1 \\citep{Friedl:2010}.', fig.pos='!h'----
plot(patterns_ts, type = "patterns")

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
log_weight = logisticWeight(alpha = -0.1, beta = 100)
matches = 
  twdtwApply(x = ts, y = patterns_ts, weight.fun = log_weight, keep=TRUE)
slotNames(matches)
show(matches)

## ----logist-time-weight, echo = FALSE, eval = TRUE, out.width=paste0(page_width/2,'in'), fig.align='center', fig.cap='Logistic time-weight function \\code{logisticWeight} with steepness \\code{alpha=-0.1} and midpoint \\code{beta=100}. The $x$ axis shows the absolute difference between two dates in days and the $y$ axis shows the time-weight \\citep{Maus:2016}.', fig.pos='!h'----
# Maximum time difference in days 
max_diff = 366/2
# Set parameters 
alpha = -0.1
beta = 100
a = 1/max_diff
# Define the logistic weight 
log_fun = logisticWeight(alpha, beta)
# Define the linear weight 
lin_fun = linearWeight(a)
# Build data.frame with linear and logistic time-weight 
Difference = 0:max_diff
df_weight = data.frame(Difference, Logistic=log_fun(Difference), Linear=lin_fun(Difference))
names(df_weight) = c("Difference", paste0("Logistic weight, alpha: ",alpha," and beta: ",beta), paste0("Linear weight, slop: ",round(a,3)))
# Reshape and plot weight curves 
df_weight = melt(df_weight, id.vars = "Difference")
names(df_weight)[-1] = c("Functions","Weight")
ggplot(df_weight, aes_string(x="Difference", y="Weight", group="Functions", linetype="Functions")) + 
  geom_line() + xlab("Time difference (days)")  + 
  theme(text = element_text(size = 10, family="Helvetica"), 
        plot.title = element_text(size = 10, family="Helvetica", face="bold"),
        axis.title = element_text(size = 10, family="Helvetica"),
        legend.position = c(.3,.85), legend.background = element_rect(fill="transparent")) +
  scale_linetype(guide_legend(title = ""))

## ----twdtw-matches, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_height/3.5, fig.align='center', fig.cap=c('The four best matches of the "soybean" pattern in the time series using a logistic time-weight. The solid black line is the long-term time series; the coloured lines are the temporal patterns; and the grey dashed lines are the respective matching points.'), fig.pos='!h'----
plot(matches, type="matches", patterns.labels="Soybean", k=4)

## ----alignments-all-patterns, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_height/2.5, fig.align='center', fig.cap=c('Alignments and dissimilarity measures of the patterns "soybean", "cotton", and "maize" to the subintervals of the long-term time series using a logistic time-weight. The solid black line is the EVI time series, and the coloured lines are the alignments of the patterns that have dissimilarity measure lower than three.'), fig.pos='!h'----
plot(matches, type="alignments", attr = "evi", threshold = 3.0)

## ----time-series-classification, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_height/2.5, fig.align='center', fig.cap=c('Classification of each 6 months periods of the time series using results of the TWDTW analysis with logistic time-weight. The solid lines are the attributes of the time series, the background colours indicate the classification of the periods.'), fig.pos='!h'----
ts_classification = twdtwClassify(x = matches, 
  from = as.Date("2009-09-01"), to = as.Date("2013-09-01"), 
  by = "6 month", overlap = 0.5)
plot(ts_classification, type="classification")

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
data_folder = system.file("lucc_MT/data", package = "dtwSat")
dir(data_folder)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
blue = brick(paste(data_folder,"blue.tif", sep = "/"))
red  = brick(paste(data_folder,"red.tif", sep = "/"))
nir  = brick(paste(data_folder,"nir.tif", sep = "/"))
mir  = brick(paste(data_folder,"mir.tif", sep = "/"))
evi  = brick(paste(data_folder,"evi.tif", sep = "/"))
ndvi = brick(paste(data_folder,"ndvi.tif", sep = "/"))
day_of_year  = brick(paste(data_folder,"doy.tif", sep = "/"))
dates = scan(paste(data_folder,"timeline", sep = "/"), what = "dates")

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
field_samples = read.csv(paste(data_folder,"samples.csv", sep = "/"))
head(field_samples, 2)
table(field_samples[["label"]])
proj_str = scan(paste(data_folder,"samples_projection", sep = "/"), 
  what = "character")
proj_str

## ---- echo = TRUE, eval = TRUE-------------------------------------------
raster_timeseries = twdtwRaster(blue, red, nir, mir, evi, ndvi, 
  timeline = dates, doy = day_of_year)

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
field_samples_ts = getTimeSeries(raster_timeseries, 
  y = field_samples, proj4string = proj_str)
field_samples_ts

## ---- echo = TRUE, eval = TRUE-------------------------------------------
temporal_patterns = 
  createPatterns(field_samples_ts, freq = 8, formula = y ~ s(x))

## ----temporal-patterns, echo = TRUE, eval = TRUE, fig.width=page_width, fig.height=page_width/1.5, fig.align='center', fig.pos='!h', fig.cap='Temporal patterns of forest, cotton-fallow, soybean-cotton, soybean-maize, and soybean-millet based on the ground truth samples.'----
plot(temporal_patterns, type = "patterns") + 
  theme(legend.position = c(.8,.25))

## ---- echo = TRUE, eval = TRUE, results = 'markup'-----------------------
log_fun = logisticWeight(alpha=-0.1, beta=50) 
twdtw_dist = twdtwApply(x = raster_timeseries, y = temporal_patterns, 
  overlap = 0.5, weight.fun = log_fun, overwrite=TRUE, format="GTiff")

## ----plot-dissmilarity2013, echo = TRUE, eval = TRUE, fig.width=page_width, fig.align='center', fig.cap='Illustration of the TWDTW dissimilarity from each temporal pattern in 2013. The blue areas are more similar to the pattern and the red areas are less similar to the pattern.', fig.pos='!h'----
plot(x = twdtw_dist, type="distance", time.levels = 6)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
land_use_maps = twdtwClassify(twdtw_dist, format="GTiff", overwrite=TRUE)

## ----plot-map, echo = TRUE, eval = TRUE, fig.width=page_width, fig.align='center', fig.cap='Land use maps for each year from 2008 to 2013.', fig.pos='!h'----
plot(x = land_use_maps, type = "maps")

## ----plot-area, echo = TRUE, eval = TRUE, fig.width=page_width, fig.align='center', fig.cap='Percentage of area for each land use class from 2008 to 2013.', fig.pos='!h'----
plot(x = land_use_maps, type = "area")

## ----plot-change, echo = TRUE, eval = TRUE, fig.width=page_width, fig.align='center', fig.cap='Gains and losses in area from the other classes. The $y$ axis shows the actual class; the positive direction of $x$ axis shows the gains and the negative direction of $x$ axis shows the losses of the classes indicated in $y$. The colors indicate from/to which classes the gains/losses belong.', fig.pos='!h'----
plot(x = land_use_maps, type = "changes")

## ----plot-dissmilarity, echo = TRUE, eval = TRUE, fig.width=page_width, fig.align='center', fig.cap='TWDTW dissimilarity measure for each pixel over each classified period. The blue areas have high confidence and the red areas have low confidence in the classification.', fig.pos='!h'----
plot(x = land_use_maps, type="distance")

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  set.seed(1)
#  partitions = splitDataset(field_samples_ts, p=0.1, times=100,
#    freq = 8, formula = y ~ s(x, bs="cc"))

## ---- echo = TRUE, eval = FALSE, results = 'markup'----------------------
#  log_fun = logisticWeight(alpha=-0.1, beta=50)
#  twdtw_res = lapply(partitions, function(x){
#    res = twdtwApply(x = x$ts, y = x$patterns, weight.fun = log_fun, n=1)
#    twdtwClassify(x = res)
#  })
#  assessment = twdtwAssess(twdtw_res)
#  head(assessment, 5)

## ---- echo = FALSE, eval = TRUE------------------------------------------
load(system.file("lucc_MT/assessment.RData", package = "dtwSat"))

## ----plot-accuracy, echo = FALSE, eval = TRUE, fig.width=page_width, fig.height=page_width/2, fig.align='center', fig.cap='User\'s Accuracy (UA) and Producer\'s Accuracy (PA) of the TWDTW method for land cover classification. The plot shows the averages and their confidence interval for 99\\%.', fig.pos='!h'----
df = melt(assessment[,-1], id="label")
df$variable = factor(df$variable, levels = c("UA", "PA"), labels = c("User's Accuracy", "Producer's Accuracy"))
ggplot(df, aes(x=label, y=value)) +
  stat_summary(fun.data="mean_cl_boot", fun.args=list(conf.int = .99), 
               width=0.5, geom="crossbar", size=0.1, fill = "gray") + 
  geom_point(size=0.2) +  facet_grid(. ~ variable) + 
  scale_y_continuous(limits = c(0,1), labels = percent, breaks = seq(0,1,.2)) + 
  xlab("") + ylab("Accuracy") + coord_flip()

## ---- echo = FALSE, eval = TRUE, results = 'asis'------------------------
assess_mean = aggregate(assessment[, c("UA","PA")], list(assessment$label), mean)
assess_sd = aggregate(assessment[, c("UA","PA")], list(assessment$label), sd)
l_names = levels(assessment$label)
names(l_names) = l_names
ic_ua = t(sapply(l_names, function(i) 100*mean_cl_boot(x = assessment$UA[assessment$label==i], conf.int = .99)))[,-1]
ic_pa = t(sapply(l_names, function(i) 100*mean_cl_boot(x = assessment$PA[assessment$label==i], conf.int = .99)))[,-1]

assess_table = data.frame(
  Class = assess_mean$Group.1,
  
  MUA = sprintf("%.2f", round(100*assess_mean$UA,2)), 
  SDUA = sprintf("(%.2f)", round(100*assess_sd$UA,2)), 
  CIUA = sprintf("[%.2f-%.2f]", round(as.numeric(ic_ua[,1]),2), round(as.numeric(ic_ua[,2]),2)), 

  MPA = sprintf("%.2f", round(100*assess_mean$PA,2)),
  SDPA = sprintf("(%.2f)", round(100*assess_sd$PA,2)),
  CIPA = sprintf("[%.2f-%.2f]", round(as.numeric(ic_pa[,1]),2), round(as.numeric(ic_pa[,2]),2))
  )

x_assess = xtable::xtable(assess_table, 
          format = tab_format, digits = 2, label = "tab:assessment", alig=c("l","c","c","c","c","c","c","c"),
          caption="User\'s and Producer\'s Accuracy of the land use classification based on TWDTW analysis. $\\mu$ is the average accuracy, $\\sigma$ the standard deviation, and CI is the confidence interval of 99\\% using 100 resampling-with-replacement.")

addtorow = list()
addtorow$pos = list(0)
addtorow$command = paste("Class & \\multicolumn{3}{c}{User's Accuracy (UA) \\%} & \\multicolumn{3}{c}{Producer's Accuracy (PA)\\%}\\\\", paste(c("","$\\mu$","$\\sigma$","CI","$\\mu$","$\\sigma$","CI"), collapse="&"), "\\\\", collapse = "")

xtable::print.xtable(x_assess, add.to.row=addtorow, include.colnames = FALSE, include.rownames = FALSE, 
             comment = FALSE, caption.placement = "bottom")

