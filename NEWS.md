# dtwSat v0.2.8

* Adds faster implementation of TWDTW for logistic weight function

* Fixes small bugs 

* Adds vignettes 

# dtwSat v0.2.7

* Adds support to user defined TWDTW weight function

* Drop support to parallel processing

* Adds a minimalist function called twdtwReduceTime that is 3x faster than twdtwApply. This function can be used for high level parallel processing implemented by users

# dtwSat v0.2.6

* Fixes error in to - from : non-numeric argument to binary operator in "twdtwAssess"

* Fixes bug in .twdtw fundtion 

* Adds function for fast time series classification "twdtw_reduce_time" (~3x faster than twdtwApply)

# dtwSat v0.2.5

* Adds dtwSat paper published on Journal of Statistical Software 

* Fixing bugs 
  
  Fix error in plotAccuracy 
  
  Generalizes twdwAssess to cases with only one map 
  
  Fixes error in getTimeSeries due to time series with only one no observation 

# dtwSat v0.2.4

* New features 

  Include the function twdtwApplyParallel for TWDTW parallel processing using the package snow 
  
  Include writeRaster for twdtwRaster class 
  
  Improve tests and documentation 
  
  Improve memory usage of twdtwApply 
  
  Improve memory usage and speed of twdtwClassify 
  
  Auto recognition of the argument "doy" to avoid naming the argument "doy = doy" 
  
* Fixing bugs
  
  Fix bug in twdtwAssess for class twdtwMatches 
  
  Fix bug in twdtwRaster 
  
# dtwSat v0.2.3

* New features 

  Register TWDTW as a distance function into package proxy 
  
* Fixing bugs

  Fix typos in plot labels 

# dtwSat v0.2.2 Release Notes

* New features

  New accuracy metrics (twdtwAssess) for classified map, including User's and Producer's accuracy, and area uncertainty. 

  Include methods for accuracy visualization (plot and LaTeX tables) 

* Update data set names 

  Rename the data sets in ordes to avoid future overwriting of functions and data sets. "example\_ts" replaced with "MOD13Q1.ts". Tthe data sets are now called:
    
    MOD13Q1.MT.yearly.patterns	Data: patterns time series
    MOD13Q1.patterns.list	Data: patterns time series
    MOD13Q1.ts	Data: An example of satellite time series
    MOD13Q1.ts.labels	Data: Labels of the satellite time series in MOD13Q1.ts
    MOD13Q1.ts.list

* Fixing bugs

  Fix bug in twdtwApply wrong sign in 'by' argument
  
  Fix bug in time index for twdtwApply-twdtwRaster
  
# dtwSat v0.2.1 Release Notes

* Fix Solaris installation errors. 

# dtwSat v0.2.0 Release Notes

* Include Fortran optimization 
 
   This version includes functions written in Fortran. 

* Obsolete features 

  The S4 class 'twdtw' no longer exists. 

* New features
 
  New S4 classes: twdtwTimeSeries, twdtwMatches, and twdtwRaster.
 
  plot methods for twdtwRaster object: 'maps', 'area', 'changes', and 'distance'.
 
  plot methods for twdtwTimeSeries objects: ''patterns'' and ''timeseries''.
 
  plot methods for twdtwMatches objects: ''paths'', ''matches'', ''alignments'', ''classification'', ''cost'', ''patterns'', and ''timeseries''.
 
  createPattern function to create temporal patterns based on set of time series.
 
  getTimeSeries extract time series from raster objects.
 
  twdtwApply apply the TWDTW analysis for raster and time series objects.


# dtwSat v0.1.1 Release Notes

* New features
 
  'normalizeQuery' new normalization feature for TWDTW
 
  'template.list' new dataset. List of template time series  
 
  arguments 'from' and 'to' in 'classifyIntervals' updated to include 'character' or 'Dates' in in the format 'yyyy-mm-dd'
 
    Align query and template by name if names not null in 'twdtw' function

* deprecated features
 
    argument 'x' from function 'waveletSmoothing' is deprecated and is scheduled to be removed in the next version. Please use 'timeseries' instead.
 
    argument 'template' from functions 'twdtw' and 'mtwdtw' is deprecated and is scheduled to be removed in the next version. Please use 'timeseries' instead.
 
  argument 'normalized' is deprecated and is scheduled to be removed in the next version from all methods 
 
  'createTimeSequence' is deprecated. Use 'getModisTimeSequence' instead.
 
  Fix function name. 'classfyIntervals' is deprecated. Use 'classifyIntervals' instead.

* Fixing bugs
 
  Fix plot intervals in plotClassify
 
  replace range(x) for range(x, na.rm=TRUE) in all methods 
 
  Bug fixed in cost matrix indexing 

 
# dtwSat v0.1.0 Release Notes

* First version of dtwSat on CRAN

# dtwSat v0.0.1 Release Notes

* Earlier dtwSat development version
