###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   Assess accuracy - 2016-02-01                              #
#                                                             #
###############################################################

assessAccuracy = function(ts.list, reference, n, p){
  
  training_sample = createDataPartition(y = reference, times = n, p = p, list = TRUE)
  
  assess_list = lapply(names(training_sample), function(i){
    print(paste("Running", i))
    I = training_sample[[i]]
    
    # Split training and validation samples 
    ts_training = ts.list[I]
    ts_validation = ts.list[-I]
    reference_training = reference[I]
    reference_validation = reference[-I]

    # Group land-use classes 
    groups = as.character(unique(reference_training))
    names(groups) = groups
    J = lapply(groups, function(x) reference_training==x )
    
    # Create temporal patterns 
    samples_list = lapply(J, function(j) ts_training[j] )
    patterns_list = lapply(samples_list, createPattern, freq=8, 
                           from="2007-09-01", to="2008-09-01",
                           formula = y ~ s(time, bs="cc"))
    
    # Apply twdtw for the validation samples 
    twdtw_results = mclapply(ts_validation, mc.cores=mc.cores, FUN=twdtw, 
                             patterns = patterns_list, weight.fun = weight.fun)
    
    # Classify twdtw results
    res = do.call("rbind", mclapply(seq_along(twdtw_results), mc.cores=mc.cores, function(k){
      from =  as.Date(paste0(format(min(index(ts_validation[[k]])),"%Y"), "-09-01"))
      to =    as.Date(paste0(format(max(index(ts_validation[[k]])),"%Y"), "-09-01"))
      pred = classifyIntervals(twdtw_results[[k]], from = from, to = to, by = "12 month", overlap = overlap)
      data.frame(
        Reference = as.character(reference_validation[k]),
        Predicted = as.character(pred[["pattern"]]), 
        stringsAsFactors = FALSE
      )
    }))
    return(res)
  })
  assess_list
}
  



