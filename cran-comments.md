## Test environments
* local Fedora release 20 (Heisenbug), R 3.2.0
* win-builder 

## REVIEWS

# v0.2.0
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  6.3Mb
  sub-directories of 1Mb or more:
    tif_MT   5.8Mb

  The 'tif' files were included in order to give an example of spatiotemporal analysis using dtwSat. This example is the smallest that is still meaningful for a spatiotemporal land use analysis. This files are used in the new vignette about land use changes analysis and in the documentation examples. 

* Fortran optimization 
  
  Two functions written in Fortran were included in this version, one to perform the accumulated cost and the second to trace back the alignment paths


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
