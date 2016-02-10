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
#   R Package dtwSat - 2016-02-10                             #
#                                                             #
###############################################################


#' @title split the time series samples and create patterns 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function splits the set of samples for 
#' training and validation.
#' 
#' @param timeseries A list of \code{\link[zoo]{zoo}} objects with the time series.
#' 
#' @param ref A vector with the identification of each time series. 
#' It must have the same length as \code{x}.
#' 
#' @param times Number of partitions to create.
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{createPattern}}.
#' 
#' 
#' @docType methods
#' @return A \code{\link[base]{list}} of three attributes: the trained \code{pattrens},
#' the validation dataset \code{x}, and the reference for the validation dataset \code{ref}.
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwAssessment}}
#' 
#' @examples
#' 
#' set.seed(1)
#' field_samples = read.csv(system.file("lucc_MT/samples.csv", package="dtwSat")) 
#' proj_str = scan(file = system.file("lucc_MT/samples_projection", package="dtwSat"), 
#'                  what = "character")
#' reference = as.character(field_samples[["class"]])
#' load(system.file("lucc_MT/ts_list.RData", package="dtwSat"))
#' 
#' x = splitted_dataset = splitDataset(timeseries = ts_list, ref = reference, 
#'        times=2, p=0.1, mc.cores=1, freq=8, from="2007-09-01", to="2008-09-01", 
#'        formula = y ~ s(time, bs="cc"))
#' 
#' names(x)
#' names(x$Resample1)
#' plotPatterns(x$Resample1$patterns)
#' 
#' @export
splitDataset = function(timeseries, ref, times=1, p=0.1, mc.cores=1, ...) {
  
  partitions = createDataPartition(y = ref, times = times, p = p, list = TRUE)
  
  res = mclapply(partitions, mc.cores = mc.cores, 
                 mc.preschedule = FALSE, function(I){
                   # Split training and validation samples 
                   ts_training = timeseries[I]
                   ts_validation = timeseries[-I]
                   ref_training = ref[I]
                   ref_validation = ref[-I]
                   
                   # Group land-use classes 
                   groups = as.character(unique(ref_training))
                   names(groups) = groups
                   J = lapply(groups, function(x) ref_training==x )
                   
                   # Create temporal patterns 
                   samples_list = lapply(J, function(j) ts_training[j] )
                   patterns_list = lapply(samples_list, createPattern, ...)
                   ts_partitions = list(patterns = patterns_list, x = ts_validation, ref = ref_validation)
                   ts_partitions
                 })
  res
}

