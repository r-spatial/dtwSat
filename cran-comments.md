## Test environments
* local Fedora release 20 (Heisenbug), R 3.2.0
* win-builder 

## REVIEWS

# v0.2.0

* Include Fortran optimization 
  
  Two functions written in Fortran were included in this version, one to perform the accumulated cost and the second to trace back the alignment paths

* Extend dtwSat class for multiple alignments 

  The dtwSat class and its methods were extended to afford multiple alignments 

* Remove dependency 

  Dependency from dtw package replaced by the import of DTW step patterns from that package 

* New features

  twdtw extended to afford multiple temporal patterns in the same call
  
  twdtw method extended to receive a weight function defined by the user
  
  plotAlignments visualization of alignments and distances for several patterns in the plot 
  
  plotMatch visualization of the matching points for several alignments in the same plot 
  
  plotGroup visualization of best classe for each time interval  

* Update documentation
 
  General review in the documentation

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
