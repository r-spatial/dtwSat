## Test environments
* local Fedora release 23, R 3.2.3
* win-builder 

## REVIEWS

# v0.2.0
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is 6.7Mb
  sub-directories of 1Mb or more:
    lucc_MT  6.2Mb

  This version includes 'tif' files for an example of spatiotemporal analysis using dtwSat. The example is small but meaningful for a spatiotemporal analysis of land use changes. The 'tif' files are used in the documentation of the package and in the vignette of the package. 

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
