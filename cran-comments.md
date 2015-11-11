## Test environments
* local Fedora release 20 (Heisenbug), R 3.2.0
* win-builder 

## REVIEWS

# v0.2.0

* Include Fortran optimization 

* Remove dtw package dependency 

# v0.1.1

* Documentation update
 
  General review in the documentation.

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
