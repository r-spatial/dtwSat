# dtwSat v0.2.0 Release Notes

* Include Fortran optimization 
  
  Two functions written in Fortran were included in this version, one to perform the accumulated cost and the second to trace back the alignment paths

* Extend dtwSat class for multiple alignments 

  The dtwSat class and its methods were extended to afford multiple alignments 

* Remove dependency 

  Dependency from dtw package replaced by the import of DTW step patterns from that package 

* New features
  
  plotLUCC method for raster object. Land use and land cover change analysis 
  
  plotPatterns method to plot temporal patterns 
  
  createPattern function to create temporal patterns from several samples time series 
  
  extractSampleTimeSeries extract time series from raster objects 
  
  twdtwApply for raster and 3D arrays

  twdtw extended to afford multiple temporal patterns in the same call
  
  twdtw method extended to receive a weight function defined by the user
  
  plotAlignments visualization of alignments and distances for several patterns in the plot 
  
  plotMatch visualization of the matching points for several alignments in the same plot 
  
  plotGroup visualization of best classe for each time interval  

* Update documentation
 
  General review in the documentation


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
