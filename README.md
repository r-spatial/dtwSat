# dtwSat

<!-- badges: start -->
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/vwmaus/dtwSat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vwmaus/dtwSat/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vwmaus/dtwSat/main.svg)](https://app.codecov.io/gh/vwmaus/dtwSat)
[![CRAN](https://www.r-pkg.org/badges/version/dtwSat)](https://cran.r-project.org/package=dtwSat)
[![Downloads](https://cranlogs.r-pkg.org/badges/dtwSat?color=brightgreen)](https://www.r-pkg.org/pkg/dtwSat)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/dtwSat)](http://www.r-pkg.org/pkg/dtwSat)
<!-- badges: end -->
  

### Time-Weighted Dynamic Time Warping for satellite image time series analysis

The package *dtwSat* provides an implementation of the Time-Weighted
Dynamic Time Warping (TWDTW) method for land cover mapping using
multi-band satellite image time series (Maus et al. 2016, 2019).
*dtwSat* provides full cycle of land cover classification using image
time series, ranging from selecting temporal patterns to visualizing,
and assessing the results.

## Installing

Install either from CRAN

``` r
install.packages("dtwSat")
```

or install the development versions from GitHub

``` r
library(devtools)
devtools::install_github("vwmaus/dtwSat")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Maus:2019" class="csl-entry">

Maus, Victor, Gilberto Camara, Marius Appel, and Edzer Pebesma. 2019.
“<span class="nocase">dtwSat</span>: Time-Weighted Dynamic Time Warping
for Satellite Image Time Series Analysis in R.” *Journal of Statistical
Software* 88 (5): 1–31. <https://doi.org/10.18637/jss.v088.i05>.

</div>

<div id="ref-Maus:2016" class="csl-entry">

Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, Fernando
M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted Dynamic
Time Warping Method for Land-Use and Land-Cover Mapping.” *IEEE Journal
of Selected Topics in Applied Earth Observations and Remote Sensing* 9
(8): 3729–39. <https://doi.org/10.1109/JSTARS.2016.2517118>.

</div>

</div>
