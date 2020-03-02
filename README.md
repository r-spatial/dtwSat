<!-- Set global env -->
<!-- 
    rmarkdown::render("README.Rmd") 
-->

dtwSat
======

[![Build
Status](https://travis-ci.org/vwmaus/dtwSat.png?branch=master)](https://travis-ci.org/vwmaus/dtwSat)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/dtwSat)](http://cran.r-project.org/package=dtwSat)
[![month](http://cranlogs.r-pkg.org/badges/dtwSat)](http://www.r-pkg.org/pkg/dtwSat)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/dtwSat)](http://www.r-pkg.org/pkg/dtwSat)

### Time-Weighted Dynamic Time Warping for satellite image time series analysis

The package *dtwSat* provides an implementation of the Time-Weighted
Dynamic Time Warping (TWDTW) method for land cover mapping using
multi-band satellite image time series (Maus et al. 2016, 2019).
*dtwSat* provides full cycle of land cover classification using image
time series, ranging from selecting temporal patterns to visualising,
and assessing the results. Bellow we show a quick demo of the package
usage.

### Install

The GitHub version requires the package *devtools*

``` r
install.packages("devtools")
devtools::install_github("vwmaus/dtwSat")
```

### Quick demo

In this quick demo we will perform a TWDTW analysis for a single time
series. The data for the analysis are a set of temporal patterns in
`MOD13Q1.patterns.list` and an example of time series in `MOD13Q1.ts` in
the Brazilian state of Mato Grosso. These time series are in `zoo`
format and come with the package installation. Suppose that we want to
know the crop type of each subinterval in following time series:

``` r
library(dtwSat)
# Create and plot object time series 
ts <- twdtwTimeSeries(MOD13Q1.ts)
class(ts)
plot(ts, type = "timeseries")
```

<img src="figure/plot-MOD13Q1.ts-ts-1.png" alt="Fig. 1. Example time series which we want to classify."  />
<p class="caption">
Fig. 1. Example time series which we want to classify.
</p>

For this region in Brazil we have a set of well known temporal patterns
derived from field observations, such that:

``` r
# Create and plot object time series 
patt <- twdtwTimeSeries(MOD13Q1.MT.yearly.patterns)
class(patt)
plot(patt, type = "patterns") 
```

<img src="figure/plot-patterns-1.png" alt="Fig. 2. Typical temporal patterns of *soybean*, *cotton*, and *maize*."  />
<p class="caption">
Fig. 2. Typical temporal patterns of *soybean*, *cotton*, and *maize*.
</p>

Using these temporal patterns we run the TWDTW analysis, such that

``` r
# Define logistic time-weight, see Maus et al. (2016)
log_fun <- logisticWeight(alpha = -0.1, beta = 50) 
# Run TWDTW analysis 
matches <- twdtwApply(x = ts, y = patt, weight.fun = log_fun, keep = TRUE) 
```

The result is a `twdtwMatches` object with all possible matches of the
patterns to the time series

``` r
class(matches)
```

    ## [1] "twdtwMatches"
    ## attr(,"package")
    ## [1] "dtwSat"

``` r
show(matches)
```

    ## An object of class "twdtwMatches"
    ## Number of time series: 1 
    ## Number of alignments: 47 
    ## Patterns labels: Cotton-fallow Forest Low vegetation Pasture Soybean-cotton Soybean-fallow Soybean-maize Soybean-millet Soybean-sunflower Water Wetland

We can use several plot methods to visualize the results of the analysis
in the `twdtwMatches` object, for example, to plot the alignments

``` r
plot(x = matches, type = "alignments", threshold = 2)
```

<img src="figure/plot-alignment-1.png" alt="Fig. 3. TWDTW alignments over time and cost (distance) in y-axis."  />
<p class="caption">
Fig. 3. TWDTW alignments over time and cost (distance) in y-axis.
</p>

to plot matching point

``` r
plot(x = matches, type = "matches", attr = "evi", patterns.labels = "Soybean-cotton", k <- 1) 
```

<img src="figure/plot-match-1.png" alt="Fig. 1. The four best matches of *soybean*."  />
<p class="caption">
Fig. 1. The four best matches of *soybean*.
</p>

to plot minimum cost paths

``` r
plot(x = matches, type = "paths", patterns.labels = "Soybean-cotton") 
```

<img src="figure/plot-path-1.png" alt="Fig. 2. The minimum cost path of the TWDTW alignment for each crop type."  />
<p class="caption">
Fig. 2. The minimum cost path of the TWDTW alignment for each crop type.
</p>

and, finally to classify the subintervals of the time series. The plot
will select the best match for each period of 6 months, i.e. the class
for each period.

``` r
plot(x = matches, type = "classification",
     from = "2009-09-01", to = "2014-08-31", 
     by = "12 month", overlap = 0.5) 
```

<img src="figure/plot-group-1.png" alt="Fig. 3. Classification using the best match for each subinterval."  />
<p class="caption">
Fig. 3. Classification using the best match for each subinterval.
</p>

### Raster time series classification

The next example shows how to classify a raster time series, i.e. the
same as we did in the quick demo but now for each pixel location. For
that we use a set of MODIS (MOD13Q1 product) images from 2007 to 2013
for a region in the Brazilian Amazon. These data is included in the
package installation. Load raster time series:

``` r
evi  <- brick(system.file("lucc_MT/data/evi.tif",  package = "dtwSat"))
ndvi <- brick(system.file("lucc_MT/data/ndvi.tif", package = "dtwSat"))
red  <- brick(system.file("lucc_MT/data/red.tif",  package = "dtwSat"))
blue <- brick(system.file("lucc_MT/data/blue.tif", package = "dtwSat"))
nir  <- brick(system.file("lucc_MT/data/nir.tif",  package = "dtwSat"))
mir  <- brick(system.file("lucc_MT/data/mir.tif",  package = "dtwSat"))
doy  <- brick(system.file("lucc_MT/data/doy.tif",  package = "dtwSat"))
```

Load the dates of the MODIS images:

``` r
timeline <- scan(system.file("lucc_MT/data/timeline", package = "dtwSat"), what = "date")
```

Build raster time series:

``` r
rts <- twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
```

Load the set of ground truth samples and projection information:

``` r
field_samples <- read.csv(system.file("lucc_MT/data/samples.csv", package = "dtwSat"))
proj_str <- scan(system.file("lucc_MT/data/samples_projection", package = "dtwSat"), what = "character")
```

We use the package *caret* to split the samples into training (10%) and
validation (90%)

``` r
library(caret)
set.seed(1)
I <- unlist(createDataPartition(field_samples$label, p = 0.1))
training_samples <- field_samples[I, ]
validation_samples <- field_samples[-I, ]
```

Extract training time series from raster time series

``` r
training_ts <- getTimeSeries(rts, y = training_samples, proj4string = proj_str)
validation_ts <- getTimeSeries(rts, y = validation_samples, proj4string = proj_str)
```

Create temporal patterns using training samples

``` r
temporal_patterns <- createPatterns(training_ts, freq = 8, formula = y ~ s(x))
```

``` r
plot(temporal_patterns, type = "patterns") 
```

<img src="figure/plot-patterns-map-1.png" alt="Fig. 4. Typical temporal patterns of *Cotton-fallow*, *Forest*, *Soybean-cotton*, *Soybean-maize*, and *Soybean-millet*."  />
<p class="caption">
Fig. 4. Typical temporal patterns of *Cotton-fallow*, *Forest*,
*Soybean-cotton*, *Soybean-maize*, and *Soybean-millet*.
</p>

Apply TWDTW analysis:

``` r
# Define logistic time-weight, see Maus et al. (2016)
log_fun <- logisticWeight(-0.1, 50)

# Run serial TWDTW analysis 
r_twdtw <- twdtwApply(x = rts, y = temporal_patterns, weight.fun = log_fun, progress = 'text')

# or Run parallel TWDTW analysis
beginCluster()
r_twdtw <- twdtwApplyParallel(x = rts, y = temporal_patterns, weight.fun = log_fun, progress = 'text')
endCluster()
```

Classify raster raster time series using the results from the TWDTW
analysis

``` r
r_lucc <- twdtwClassify(r_twdtw, progress = 'text')
```

Visualising the results.

Land cover maps

``` r
plot(x = r_lucc, type = "maps")
```

<img src="figure/plot-maps-1.png" alt="Fig. 5. Land cover maps based on TWDTW analysis."  />
<p class="caption">
Fig. 5. Land cover maps based on TWDTW analysis.
</p>

Land cover area for each class over time

``` r
plot(x = r_lucc, type = "area")
```

<img src="figure/plot-area-1.png" alt="Fig. 6. Land cover area based on TWDTW analysis."  />
<p class="caption">
Fig. 6. Land cover area based on TWDTW analysis.
</p>

Land cover changes over time (gains and losses from/to classes)

``` r
plot(x = r_lucc, type = "changes")
```

<img src="figure/plot-changes-1.png" alt="Fig. 7. Land cover changes based on TWDTW analysis."  />
<p class="caption">
Fig. 7. Land cover changes based on TWDTW analysis.
</p>

We use the validation samples to compute the metrics for accuracy
assessment.

``` r
twdtw_assess <- twdtwAssess(object = r_lucc, y = validation_samples, 
  proj4string = proj_str, conf.int = .95) 
show(twdtw_assess)
```

    ## An object of class "twdtwAssessment"
    ## Number of classification intervals: 6 
    ## Accuracy metrics summary
    ## 
    ## Overall
    ## Accuracy      Var       sd      ci* 
    ##  0.95831  0.00011  0.01037  0.02033 
    ## 
    ## User's
    ##                Accuracy     Var    sd   ci*
    ## Cotton-fallow      0.95 0.00071 0.027 0.052
    ## Forest             1.00 0.00000 0.000 0.000
    ## Soybean-cotton     1.00 0.00000 0.000 0.000
    ## Soybean-maize      0.91 0.00063 0.025 0.049
    ## Soybean-millet     1.00 0.00000 0.000 0.000
    ## unclassified       1.00 0.00000 0.000 0.000
    ## 
    ## Producer's
    ##                Accuracy     Var    sd  ci*
    ## Cotton-fallow      1.00 0.00000 0.000 0.00
    ## Forest             1.00 0.00000 0.000 0.00
    ## Soybean-cotton     0.67 0.00524 0.072 0.14
    ## Soybean-maize      1.00 0.00000 0.000 0.00
    ## Soybean-millet     0.92 0.00094 0.031 0.06
    ## unclassified       1.00 0.00000 0.000 0.00
    ## 
    ## Area and uncertainty
    ##                 Mapped Adjusted     ci*
    ## Cotton-fallow  4.8e+07  4.6e+07 2518092
    ## Forest         7.4e+07  7.4e+07       0
    ## Soybean-cotton 1.6e+07  2.4e+07 5047004
    ## Soybean-maize  1.2e+08  1.1e+08 6036628
    ## Soybean-millet 6.1e+07  6.7e+07 4373953
    ## unclassified   0.0e+00  0.0e+00       0
    ## 
    ## * 95 % confidence interval

Visualizing User’s and Producer’s accuracy

``` r
plot(twdtw_assess, type = "accuracy")
```

<img src="figure/plot-users-prodcucers-1.png" alt="Fig. 8. User's and Producer's accuracy."  />
<p class="caption">
Fig. 8. User’s and Producer’s accuracy.
</p>

Visualizing area uncertainty

``` r
plot(twdtw_assess, type = "area")
```

<img src="figure/plot-area-uncertainty-1.png" alt="Fig. 9. Area uncertainty."  />
<p class="caption">
Fig. 9. Area uncertainty.
</p>

For further discussion on the package see the \[vignettes\]\[Vignettes\]
and if you want to learn more about the TWDTW method (see, Maus et al.
2016, and @Maus:2019).

### References

Maus, Victor, Gilberto Camara, Marius Appel, and Edzer Pebesma. 2019.
“dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image Time
Series Analysis in R.” *Journal of Statistical Software* 88 (5): 1–31.
<https://doi.org/10.18637/jss.v088.i05>.

Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, Fernando
M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted Dynamic
Time Warping Method for Land-Use and Land-Cover Mapping.” *IEEE Journal
of Selected Topics in Applied Earth Observations and Remote Sensing* 9
(8): 3729–39. <https://doi.org/10.1109/JSTARS.2016.2517118>.
