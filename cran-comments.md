## Test environments
* local Fedora 23 (64-bit), R 3.3.0
* win-builder 

## REVIEWS

# v0.2.0
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  7.3Mb
  sub-directories of 1Mb or more:
    lucc_MT   6.1Mb

  This version includes 'tif' files used to reproduce the examples in the vignette and documentation. The data set is the smallest that is still meaningful for a spatiotemporal analysis of land cover changes. This examples and data sets are crutial for the user to learn how to use the package. 

* Optimization 
 
  + Fortran code for optimization. 

# v0.1.0

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
