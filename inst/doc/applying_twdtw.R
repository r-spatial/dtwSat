## ---- echo = FALSE, eval = TRUE-------------------------------------------
options(replace.assign = TRUE, width = 76, prompt = "R> ")

## ---- echo=FALSE, cache = FALSE-------------------------------------------
library(knitr)
opts_chunk$set(
  prompt     = TRUE,
  width      = 76,
  warning = FALSE,
  message = FALSE,
  error = FALSE,
  results = "hide",
  cache.path = "./cache/",
  cache = FALSE,
  comment = ""
)
#knit_hooks$set(purl = hook_purl) # Save chunks to Rscript file  

###### Test 3
# Work around encoding issue on Windows
other_user <- Sys.getenv("USER") != "maus"
set.seed(20160708L)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass{jss}\\usepackage{siunitx}\\usepackage{latexsym}")
if (other_user) {
  opts_chunk$set(dev="png")
} else {
  opts_chunk$set(dev="pdf")
}
##### End test 3

#options(tikzDocumentDeclaration = "\\documentclass{jss}\\usepackage{siunitx}\\usepackage{latexsym}")

## ---- echo=FALSE, eval = TRUE, cache = FALSE------------------------------
library(dtwSat)
library(ggplot2)
library(scales)
library(Hmisc)
library(reshape2)

new_theme = theme_get()
new_theme$text$family = "Helvetica"
new_theme$text$size = 8
old_theme = theme_set(new_theme)

tab_format = "latex"
page_width = 5.590551#in 14.2#cm
page_height = 9.173228#in 23.3#cm

## ----twdtw-example, echo = FALSE, eval = TRUE, fig.width = page_width, fig.height = page_height/4, fig.align = 'center', fig.cap = 'Matches of the known temporal pattern to subintervals of the long-term time series. The solid black line is the long-term time series, the colored lines are the different matches of the same pattern, and the gray dashed lines are the matching points.', fig.pos = '!h'----
n    <- 4
ts   <- twdtwTimeSeries(MOD13Q1.ts, labels = "Time series")
patt <- twdtwTimeSeries(MOD13Q1.patterns.list$Cotton, labels = "Matches")
log_fun <- logisticWeight(alpha = -0.1, beta = 100)
mat     <- twdtwApply(x = ts, y = patt, weight.fun = log_fun, keep = TRUE, n = n)
df_dist <- mat[[1]]
df_dist$label <- paste("Distance:", round(df_dist$distance, 2))
df_dist$y     <- 1.8 
plotMatches(mat, attr = "evi", k = n) + 
  ylab("Time series                  Pattern") + 
  geom_text(data = df_dist, mapping = aes_string(x = 'to', y = 'y', label = 'label'),
            size = 2, family="Helvetica") + 
  theme(legend.position = "none")

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
library(dtwSat)
names(MOD13Q1.patterns.list)
head(MOD13Q1.ts, n = 2)

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
ts <- twdtwTimeSeries(MOD13Q1.ts, labels = "Time series")
patterns_ts <- twdtwTimeSeries(MOD13Q1.patterns.list)
patterns_ts

## ----example-timeseries, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_height / 3, fig.align = 'center', fig.cap = 'Example of time series based on MODIS product MOD13Q1 \\citep{Friedl:2010}. The labels of the phenological cycle are shown in the plot.', fig.pos = '!h'----
plot(ts, type = "timeseries") + annotate(geom = "text", 
  x = MOD13Q1.ts.labels$from + 90, y = 0.98, 
  label = MOD13Q1.ts.labels$label, size = 2)

## ----temporal-patterns-soy-cot-mai, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_height / 3.5, fig.align = 'center', fig.cap = 'Temporal patterns of soybean, cotton, and maize based on MODIS product MOD13Q1 \\citep{Friedl:2010}.', fig.pos = '!h'----
plot(patterns_ts, type = "patterns")

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
log_weight <- logisticWeight(alpha = -0.1, beta = 100)
matches <- 
  twdtwApply(x = ts, y = patterns_ts, weight.fun = log_weight, keep = TRUE)
slotNames(matches)
show(matches)

## ----logist-time-weight, echo = FALSE, eval = TRUE, out.width = paste0(page_width / 2, 'in'), fig.align = 'center', fig.cap = 'Logistic time-weight function \\code{logisticWeight} with steepness \\code{alpha = -0.1} and midpoint \\code{beta = 100}. The $x$ axis shows the absolute difference between two dates in days and the $y$ axis shows the time-weight \\citep{Maus:2016}.', fig.pos = '!h'----
# Maximum time difference in days 
max_diff <- 366/2
# Set parameters 
alpha <- -0.1
beta  <- 100
a <- 1 / max_diff
# Define the logistic weight 
log_fun <- logisticWeight(alpha, beta)
# Define the linear weight 
lin_fun <- linearWeight(a)
# Build data.frame with linear and logistic time-weight 
Difference <- 0:max_diff
df_weight <- data.frame(Difference, Logistic = log_fun(Difference), Linear = lin_fun(Difference))
names(df_weight) <- c("Difference", paste0("Logistic weight, alpha: ", alpha, " and beta: ", beta),
                      paste0("Linear weight, slop: ", round(a, 3)))
# Reshape and plot weight curves 
df_weight <- melt(df_weight, id.vars = "Difference")
names(df_weight)[-1] <- c("Functions", "Weight")
ggplot(df_weight, aes_string(x="Difference", y="Weight", group="Functions", linetype="Functions")) + 
  geom_line() + xlab("Time difference (days)")  + 
  theme(text = element_text(size = 10, family="Helvetica"), 
        plot.title = element_text(size = 10, family="Helvetica", face="bold"),
        axis.title = element_text(size = 10, family="Helvetica"),
        legend.position = c(.3,.85), legend.background = element_rect(fill="transparent")) +
  scale_linetype(guide_legend(title = ""))

## ----twdtw-matches, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_height / 3.5, fig.align = 'center', fig.cap = c('The four best matches of the "soybean" pattern in the time series using a logistic time-weight. The solid black line is the long-term time series; the coloured lines are the temporal patterns; and the grey dashed lines are the respective matching points.'), fig.pos = '!h'----
plot(matches, type = "matches", patterns.labels = "Soybean", k = 4)

## ----alignments-all-patterns, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_height / 2.5, fig.align = 'center', fig.cap = c('Alignments and dissimilarity measures of the patterns "soybean", "cotton", and "maize" to the subintervals of the long-term time series using a logistic time-weight. The solid black line is the EVI time series, and the coloured lines are the alignments of the patterns that have dissimilarity measure lower than three.'), fig.pos = '!h'----
plot(matches, type = "alignments", attr = "evi", threshold = 3.0)

## ----time-series-classification, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_height / 2.8, fig.align = 'center', fig.cap = c('Classification of each 6 months periods of the time series using results of the TWDTW analysis with logistic time-weight. The solid lines are the attributes of the time series, the background colours indicate the classification of the periods.'), fig.pos = '!ht'----
ts_classification <- twdtwClassify(x = matches, 
  from = as.Date("2009-09-01"), to = as.Date("2013-09-01"), 
  by = "6 month", overlap = 0.5)
plot(ts_classification, type = "classification") 

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
data_folder <- system.file("lucc_MT/data", package = "dtwSat")
dir(data_folder)

## ---- echo = TRUE, eval = TRUE--------------------------------------------
blue <- brick(paste(data_folder, "blue.tif", sep = "/"))
red  <- brick(paste(data_folder,  "red.tif", sep = "/"))
nir  <- brick(paste(data_folder,  "nir.tif", sep = "/"))
mir  <- brick(paste(data_folder,  "mir.tif", sep = "/"))
evi  <- brick(paste(data_folder,  "evi.tif", sep = "/"))
ndvi <- brick(paste(data_folder, "ndvi.tif", sep = "/"))
day_of_year <- brick(paste(data_folder, "doy.tif", sep = "/"))
dates <- scan(paste(data_folder, "timeline", sep = "/"), what = "dates")

## ---- echo = TRUE, eval = TRUE--------------------------------------------
raster_timeseries <- twdtwRaster(blue, red, nir, mir, evi, ndvi, 
  timeline = dates, doy = day_of_year)

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
field_samples <- read.csv(paste(data_folder, "samples.csv", sep = "/"))
head(field_samples, 5)
table(field_samples[["label"]])
proj_str <- scan(paste(data_folder, "samples_projection", sep = "/"), 
  what = "character")
proj_str

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
field_samples_ts <- getTimeSeries(raster_timeseries, 
  y = field_samples, proj4string = proj_str)
field_samples_ts

## ---- echo = TRUE, eval = FALSE, results = 'markup'-----------------------
#  set.seed(1)
#  log_fun <- logisticWeight(alpha = -0.1, beta = 50)
#  cross_validation <- twdtwCrossValidate(field_samples_ts,
#    times = 100, p = 0.1, freq = 8, formula = y ~ s(x, bs = "cc"),
#    weight.fun = log_fun)

## ---- echo = FALSE, eval = TRUE-------------------------------------------
load(system.file("lucc_MT/cross_validation.RData", package = "dtwSat"))

## ---- echo = FALSE, eval = TRUE, results = 'asis'-------------------------
twdtwXtable(cross_validation, conf.int = .95, digits = 2, caption = "\\label{tab:cross-validation} User\'s and producer\'s accuracy of the TWDTW cross-validation using 100 resampling-with-replacement. The table shows the standard deviation ($\\sigma$) and the 95\\% confidence interval (ci) of the mean ($\\mu$).", comment = FALSE, caption.placement = "bottom", table.placement="!ht", show.footnote = FALSE)

## ----plot-accuracy, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_width / 2, fig.align = 'center', fig.cap = 'User\'s and producer\'s accuracy of the TWDTW cross-validation using 100 resampling-with-replacement. The plot shows the 95\\% confidence interval of the mean.', fig.pos = '!ht'----
plot(cross_validation, conf.int = .95)

## ---- echo = TRUE, eval = TRUE--------------------------------------------
library(caret) 
set.seed(1)
I <- unlist(createDataPartition(field_samples[ , "label"], p = 0.1))
training_ts <- subset(field_samples_ts, I)
validation_samples <- field_samples[-I, ]

## ---- echo = TRUE, eval = TRUE--------------------------------------------
temporal_patterns <- 
  createPatterns(training_ts, freq = 8, formula = y ~ s(x))

## ----temporal-patterns, echo = TRUE, eval = TRUE, fig.width = page_width, fig.height = page_width / 1.3, fig.align = 'center', fig.pos = '!h', fig.cap = 'Temporal patterns of Forest, Cotton-fallow, Soybean-cotton, Soybean-maize, and Soybean-millet based on the ground truth samples.'----
plot(temporal_patterns, 
  type = "patterns") + theme(legend.position = c(.8, .21))

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
log_fun <- logisticWeight(alpha = -0.1, beta = 50) 
twdtw_dist <- twdtwApply(x = raster_timeseries, y = temporal_patterns, 
  overlap = 0.5, weight.fun = log_fun, overwrite = TRUE, format = "GTiff")

## ----plot-dissmilarity2008, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'Illustration of the TWDTW dissimilarity from each temporal pattern in 2008. The blue areas are more similar to the pattern and the red areas are less similar to the pattern.', fig.pos = '!ht'----
plot(x = twdtw_dist, type = "distance")

## ---- echo = TRUE, eval = TRUE, results = 'markup'------------------------
land_cover_maps <- 
  twdtwClassify(twdtw_dist, format = "GTiff", overwrite = TRUE)

## ----plot-map, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'Land cover maps for each year from 2008 to 2013.', fig.pos = '!h'----
plot(x = land_cover_maps, type = "maps")

## ----plot-area, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'Percentage of area for each land cover class from 2008 to 2013.', fig.pos = '!h'----
plot(x = land_cover_maps, type = "area")

## ----plot-change, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'Gains and losses in area from the other classes. The $y$ axis shows the actual class; the positive direction of $x$ axis shows the gains and the negative direction of $x$ axis shows the losses of the classes indicated in $y$. Colors indicate which class the gain/losse belong.', fig.pos = '!h'----
plot(x = land_cover_maps, type = "changes")

## ----plot-dissmilarity, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'TWDTW dissimilarity measure for each pixel over each classified period. The blue areas have high confidence and the red areas have low confidence in the classification.', fig.pos = '!h'----
plot(x = land_cover_maps, type = "distance")

## ---- echo = TRUE, eval = TRUE--------------------------------------------
maps_assessment <- twdtwAssess(land_cover_maps, y = validation_samples, 
  proj4string = proj_str, conf.int = .95)

## ---- echo = FALSE, eval = TRUE, results = 'asis'-------------------------
twdtwXtable(maps_assessment, table.type = "errormatrix", digits = 0, rotate.col = TRUE, caption = "\\label{tab:map-error-matrix}Error matrix of the map classification based on TWDTW analysis. The area is in the map unit, in this case $m^2$. $w$ is the proportion of area mapped for each class.", comment = FALSE, caption.placement = "bottom", table.placement = "!ht", show.footnote = FALSE)

## ---- echo = FALSE, eval = TRUE, results = 'asis'-------------------------
twdtwXtable(maps_assessment, table.type = "accuracy", show.prop = TRUE, digits = 2, rotate.col = TRUE, caption = "\\label{tab:map-accuracy}Accuracy and error matrix in proportion of area of the classified map. * 95\\% confidence interval.", comment = FALSE, caption.placement = "bottom", table.placement = "!ht", show.footnote = FALSE)

## ----plot-map-incorrect-samples, echo = TRUE, eval = TRUE, fig.width = page_width, fig.align = 'center', fig.cap = 'Incorrect classified samples.', fig.pos = '!ht'----
plot(x = maps_assessment, type = "map", samples = "incorrect")

## ---- echo = FALSE, eval = TRUE, results = 'asis'-------------------------
twdtwXtable(maps_assessment, table.type = "area", digits = 0, rotate.col = TRUE, caption = "\\label{tab:map-adjusted-area}Mapped and adjusted, accumulated over the whole period, i.e., the sum of the maps from 2008 to 2013. The area is in the map unit, in this case $m^2$. * 95\\% confidence level.", comment = FALSE, caption.placement = "bottom", table.placement = "!ht", show.footnote = FALSE)

## ----plot-area-and-uncertainty, echo = FALSE, eval = TRUE, fig.width = page_width, fig.height = page_height / 2.7, fig.align = 'center', fig.cap = 'Mapped and adjusted, accumulated over the whole period, i.e., the sum from the sum of the maps from 2008 to 2013. The area is in the map unit, in this case $m^2$.', fig.pos = '!ht'----
plot(x = maps_assessment, type = "area", perc = FALSE) + 
  ylab(expression(paste("Area", m^2, sep = " ")))

