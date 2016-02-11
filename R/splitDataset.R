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
#' training and validation. The function uses stratified sampling. 
#' Then a simple random sampling is applied within each stratum.
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
#' @param mc.cores The number of cores to use, See \code{\link[parallel]{mclapply}} 
#' for details.
#' 
#' @param p the percentage of data that goes to training. 
#' See \code{\link[caret]{createDataPartition}} for details.
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
#' \dontrun{
#' set.seed(1)
#' data_folder = system.file("lucc_MT/data", package = "dtwSat")
#' field_samples = read.csv(paste(data_folder,"samples.csv", sep = "/"))
#' proj_str = scan(paste(data_folder,"samples_projection", sep = "/"), 
#'            what = "character")
#' reference = as.character(field_samples[["class"]])
#' load(system.file("lucc_MT/ts_list.RData", package="dtwSat"))
#' 
#' x = splitted_dataset = splitDataset(timeseries = ts_list, ref = reference, 
#'     times=2, p=0.1, mc.cores=1, freq=8, from="2007-09-01", to="2008-09-01", 
#'     formula = y ~ s(time, bs="cc"))
#' 
#' names(x)
#' names(x$Resample1)
#' plotPatterns(x$Resample1$patterns)
#' }
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

