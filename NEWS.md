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
