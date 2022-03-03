getMatchingDates <- function(x){
  
  best_aligs <- x$internals$alignments[x$internals$alignments[,6]==1,]
  best_aligs <- best_aligs[order(best_aligs[,1]),]
  
  out <- lapply(1:nrow(best_aligs), function(i){
    idx <- as.data.frame(.tracepath(dm = x$internals$internals[[i]]$DM,
                                    step.matrix = x$internals$internals[[i]]$SM,
                                    jmin = best_aligs[i,3]))
    idx$patternDates <- x$internals$internals[[i]]$patternDates[idx$index1]
    idx$tsDates <- x$internals$internals[[i]]$tsDates[idx$index2]
    return(idx)
    })
  
  return(out)

}
