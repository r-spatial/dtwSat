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
    out = timeSeriesAnalysis(query, template, theta=THETA, span=1,
                             normalize=TRUE)
    return(out)
  })
  names(res) = names(TemporalPatterns.list)
  return(res)
}
