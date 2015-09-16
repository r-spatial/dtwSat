#' @title Splits the total number of processes
#' 
#' @description A recursive function to splits the total 
#' number of processes according to the threads size.
#' 
#' @param n An integer for the total number of processes. 
#' @param thread.size An integer for the thread size.
#' @docType methods
#' @export
computeThreadsSize = function(n, thread.size){
  if( n < 1.5*thread.size || thread.size < 1)
    return(n)
  return( c(thread.size, computeThreadsSize(n-thread.size, thread.size)) )
}




#' @title Computes DTW distance
#' 
#' @description The function computes the whole possible alignments 
#' for a set of patterns in one time series. 
#' 
#' @param template A zoo object with the template time series. 
#' It must be larger than the query.
#' @param TemporalPatterns.list A list of patterns. Each node of the list 
#' has two columns data.frame with the dates and the value. 
#' @docType methods
#' @export
computeDTWForAllPatterns = function(template, TemporalPatterns.list, ... ){
  res = lapply(seq_along(TemporalPatterns.list), function(j){
    patternName = names(TemporalPatterns.list)[j]
    x = as.numeric(TemporalPatterns.list[[j]][,2])
    tx = as.Date(TemporalPatterns.list[[j]][,1], origin="1970-01-01")
    query = zoo(x, tx)
    out = timeSeriesAnalysis(query, template, theta=THETA, span=SPAN,
                             normalize=NORMALIZE, satStat=SATSTAT)
    return(out)
  })
  names(res) = names(TemporalPatterns.list)
  return(res)
}
