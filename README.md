# dtwSat

<!-- badges: start -->
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/r-spatial/dtwSat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-spatial/dtwSat/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vwmaus/dtwSat/main.svg)](https://app.codecov.io/gh/vwmaus/dtwSat)
[![CRAN](https://www.r-pkg.org/badges/version/dtwSat)](https://cran.r-project.org/package=dtwSat)
[![Downloads](https://cranlogs.r-pkg.org/badges/dtwSat?color=brightgreen)](https://www.r-pkg.org/pkg/dtwSat)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/dtwSat)](http://www.r-pkg.org/pkg/dtwSat)
[![cran checks](https://badges.cranchecks.info/worst/dtwSat.svg)](https://cran.r-project.org/web/checks/check_results_dtwSat.html)
[![status](https://tinyverse.netlify.com/badge/dtwSat)](https://CRAN.R-project.org/package=dtwSat)

<!-- badges: end -->
  
Provides a robust approach to land use mapping using multi-dimensional 
(multi-band) satellite image time series. By leveraging the Time-Weighted Dynamic 
Time Warping (TWDTW) distance metric in tandem with a 1 Nearest-Neighbor (1-NN) Classifier,
this package offers functions to produce land use maps based on distinct seasonality patterns, 
commonly observed in the phenological cycles of vegetation. The approach is described in 
Maus et al. (2016) and Maus et al. (2019).
A primary advantage of TWDTW is its capability to handle irregularly sampled and noisy time series, 
while also requiring minimal training sets. The package includes tools for training the 1-NN-TWDTW model, 
visualizing temporal patterns, producing land use maps, and visualizing the results.

## Getting Started

You can install dtwSat from CRAN using the following command:

``` r
install.packages("dtwSat")
```

Alternatively, you can install the development version from GitHub:

``` r
devtools::install_github("r-spatial/dtwSat")
```

After installation, you can read the vignette for a quick start guide:

``` r
vignette("landuse-mapping", "dtwSat")
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
