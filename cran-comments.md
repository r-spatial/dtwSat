# Test environments

* win-builder 
  devtools::check_win_release()
  devtools::check_win_devel()
  devtools::check_win_oldrelease()

* R-hub 
  rhub::check_for_cran(check_args = '--as-cran')
  rhub::check_for_cran(check_args = '--as-cran', valgrind = TRUE)

* Local Ubuntu 22.04.1 LTS x86_64-pc-linux-gnu (64-bit), R 4.3.1 (2023-06-16)
  devtools::check(args = '--as-cran')
  devtools::submit_cran()
  
# REVIEWS

## v1.0.0

* Major release that removes obsolete dependencies, such as raster, rgdal, and sp.

* Substantially reduced the number of dependencies.

* Introduced sf and stars for spatial data handling.

* Provided a workflow compatible with other image classification workflows.




## v0.2.8

* Fixes errors from https://cran.r-project.org/web/checks/check_results_dtwSat.html

* Speed improvements

## v0.2.7

* Fixes error in cost TWDTW weighting function

* Drop support to parallel processing

* Drop projection method for twdtwRaster class as it was used only internally 

* Fixes errors from https://cran.r-project.org/web/checks/check_results_dtwSat.html

## v0.2.6

* Fixes warnings from https://cran.r-project.org/web/checks/check_results_dtwSat.html

## v0.2.5

* The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.

## v0.2.4

## v0.2.3

* Fix check error 
   checking re-building of vignette outputs ... [1s/1s] WARNING 
   Error in re-building vignettes: 
   ... 

## v0.2.2

## v0.2.1

* Fix Solaris installation errors. 
    Replacing the GNU extension ISNAN with pure Fortran code to check for NAN. 

## v0.2.0

* Single quote software names in the Descriptoin.
  
    - Using dtwSat the user 
    + Using 'dtwSat' the user 

* Please reduce to < 5 MB, I do not believe this is not possible. 
  
    gdal_translate COMPRESS=DEFLATE reduced the size of the tiff files. 
  
* build_win() latex compilation error: 

I have successfully built and compiled the latex of the vignette using a personal Windows machine. However, 'build_win' gives an error while compiling the latex of vignette. The error message given by the server is unclear to me:

    * checking re-building of vignette outputs ... WARNING

    Error in re-building vignettes:

    ...

    Error: processing vignette 'applying_twdtw.Rmd' failed with diagnostics:

    Failed to compile applying_twdtw.tex.

    Execution halted





There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  7.3Mb
  sub-directories of 1Mb or more:
    lucc_MT   6.1Mb

  This version includes 'tif' files used to reproduce the examples in the vignette and documentation. The data set is the smallest that is still meaningful for a spatiotemporal analysis of land cover changes. This examples and data sets are crutial for the user to learn how to use the package. 

* Optimization 
 
  + Fortran code for optimization. 

## v0.1.0

* authors / copyright holder

  - C code removed from the package 
 
* checking R code for possible problems ... NOTE
plotCostMatrix: no visible global function definition for 'gray.colors'

  + importFrom("grDevices", "gray.colors") added to NAMESPACE

* DESCRIPTION file

  - Description: The dtwSat provides ...
  + Description: Provides ...

* Functions documentation 

  General review in the documentation.


## R CMD check --as-cran results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Victor Maus <vwmaus1@gmail.com>'
  New submission

  This is my first submission.
