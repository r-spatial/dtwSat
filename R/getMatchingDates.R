getMatchingDates <- function(x){
  
  best_aligs <- x$internals$alignments[x$internals$alignments[,6]==1,,drop=FALSE]
  best_aligs <- best_aligs[order(best_aligs[,1]),,drop=FALSE]
  
  out <- lapply(1:nrow(best_aligs), function(i){
    ts_id <- best_aligs[i,5]
    idx <- as.data.frame(.tracepath(dm = x$internals$internals[[ts_id]]$DM,
                                    step.matrix = x$internals$internals[[ts_id]]$SM,
                                    jmin = best_aligs[i,3]))
    idx$patternDates <- x$internals$internals[[ts_id]]$patternDates[idx$index1]
    idx$tsDates <- x$internals$internals[[ts_id]]$tsDates[idx$index2]
    return(idx)
    })
  
  return(out)

}
