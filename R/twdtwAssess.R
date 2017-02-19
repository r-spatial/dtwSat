
setGeneric("twdtwAssess", 
           def = function(object, ...) standardGeneric("twdtwAssess")
)

#' @title Assess TWDTW classification 
#' @name twdtwAssess
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#'  
#' @description Performs an accuracy assessment 
#' of the classified maps. The function returns Overall Accuracy, 
#' User's Accuracy, Produce's Accuracy, error matrix (confusion matrix),
#' and estimated area according to [1-2]. The function returns the metrics 
#' for each time interval and a summary considering all classified intervals. 
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwRaster}} resulting from 
#' the classification, i.e. \code{\link[dtwSat]{twdtwClassify}}.
#' The argument can also receive an error matrix (confusion matrix) in using the classes 
#' \code{\link[base]{data.frame}} or \code{\link[base]{table}}. In this case the user 
#' must inform the area for each class to the argument \code{area}. 
#' 
#' @param area a numeric vector with the area for each class if the argument \code{object}
#' is an error matrix (confusion matrix). 
#'
#' @param y a \code{\link[base]{data.frame}} whose attributes are: longitude, 
#' latitude, the start ''from'' and the end ''to'' of the time interval 
#' for each sample. This can also be a \code{\link[sp]{SpatialPointsDataFrame}} 
#' whose attributes are the start ''from'' and the end ''to'' of the time interval.
#' If missing ''from'' and/or ''to'', they are set to the time range of the 
#' \code{object}. 
#' 
#' @param id.labels a numeric or character with an column name from \code{y} to 
#' be used as samples labels. Optional.
#' 
#' @param labels character vector with time series labels. For signature 
#' \code{\link[dtwSat]{twdtwRaster}} this argument can be used to set the 
#' labels for each sample in \code{y}, or it can be combined with \code{id.labels} 
#' to select samples with a specific label.
#' 
#' @param proj4string projection string, see \code{\link[sp]{CRS-class}}. Used 
#' if \code{y} is a \code{\link[base]{data.frame}}.
#' 
#' @param conf.int specifies the confidence level (0-1).
#' 
#' @references 
#' [1] Olofsson, P., Foody, G.M., Stehman, S.V., Woodcock, C.E. (2013). 
#' Making better use of accuracy data in land change studies: Estimating 
#' accuracy and area and quantifying uncertainty using stratified estimation. 
#' Remote Sensing of Environment, 129, pp.122-131.
#' 
#' @references 
#' [2] Olofsson, P., Foody G.M., Herold M., Stehman, S.V., Woodcock, C.E., Wulder, M.A. (2014)
#' Good practices for estimating area and assessing accuracy of land change. Remote Sensing of 
#' Environment, 148, pp. 42-57.
#'
#' @seealso \code{\link[dtwSat]{twdtwClassify}},  
#' \code{\link[dtwSat]{twdtwAssessment}}, and
#' \code{\link[dtwSat]{twdtwXtable}}.
#' 
#' @examples 
#' \dontrun{
#' 
#' # Create raster time series
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' 
#' # Read fiels samples 
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", 
#'                 package="dtwSat"), what = "character")
#' 
#' # Split samples for training (10%) and validation (90%) using stratified sampling 
#' library(caret) 
#' set.seed(1)
#' I = unlist(createDataPartition(field_samples$label, p = 0.1))
#' training_samples = field_samples[I,]
#' validation_samples = field_samples[-I,]
#' 
#' # Create temporal patterns 
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))
#' 
#' # Run TWDTW analysis for raster time series 
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' r_twdtw = twdtwApply(x=rts, y=temporal_patterns, weight.fun=log_fun, format="GTiff", 
#'                      overwrite=TRUE)
#'                      
#' # Classify raster based on the TWDTW analysis 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' plot(r_lucc)
#' 
#' # Assess classification 
#' twdtw_assess = twdtwAssess(object = r_lucc, y = validation_samples, 
#'                            proj4string = proj_str, conf.int=.95) 
#' twdtw_assess
#' 
#' # Plot assessment 
#' plot(twdtw_assess, type="accuracy")  
#' plot(twdtw_assess, type="area") 
#' 
#' # Create latex tables 
#' twdtwXtable(twdtw_assess, table.type="matrix")
#' twdtwXtable(twdtw_assess, table.type="accuracy")
#' twdtwXtable(twdtw_assess, table.type="area")
#' 
#' }
NULL

#' @aliases twdtwAssess-twdtwRaster
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' @export
setMethod(f = "twdtwAssess", signature = "twdtwRaster",
          definition = function(object, y, labels=NULL, id.labels=NULL, proj4string=NULL, conf.int=.95) 
            twdtwAssess.twdtwRaster(object, y, labels, id.labels, proj4string, conf.int))

#' @aliases twdtwAssess-data.frame
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' 
#' @examples 
#' 
#' # Total mapped area by class. Data from [1]
#' area = c(A = 22353, B = 1122543, C = 610228) 
#' 
#' # Error matrix, columns (Reference) rows (Map)
#' x = data.frame(
#'               rbind(
#'               c( 97,	 0,   3),
#'               c(  3, 279,  18),
#'               c(  2,   1,  97)
#'    ))
#'
#' names(x) = names(area)
#' rownames(x) = names(area)
#' 
#' table_assess = twdtwAssess(x, area, conf.int = .95)
#' 
#' table_assess
#' 
#' plot(table_assess, type="area", perc=FALSE)
#' 
#' @export
setMethod(f = "twdtwAssess", signature = "data.frame",
          definition = function(object, area, conf.int=.95) 
            twdtwAssess.table(object, area, conf.int))

#' @aliases twdtwAssess-twdtwRaster
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' @export
setMethod(f = "twdtwAssess", signature = "table",
          definition = function(object, area, conf.int=.95) 
            twdtwAssess(as.data.frame.matrix(object), area, conf.int))

twdtwAssess.table = function(object, area, conf.int){
  
  if(ncol(object)!=nrow(object))
    stop("object has have the same number of rows and columns")
  
  if(nrow(object)!=length(area))
    stop("area must have length equal to the number of rows in object")
  
  accuracy = .twdtwAssess(object, area, conf.int)
  
  new("twdtwAssessment", accuracySummary=accuracy, accuracyByPeriod=accuracy)
}
  
twdtwAssess.twdtwRaster = function(object, y, labels, id.labels, proj4string, conf.int){
  
  # Check control points 
  y = .adjustLabelID(y, labels, id.labels)
  if(!"from"%in%names(y))
    stop("samples starting date not found, the argument 'y' must have a column called 'from'")
  if(!"to"%in%names(y))
    stop("samples ending date not found, the argument 'y' must have a column called 'to'")
  y = .toSpatialPointsDataFrame(y, object, proj4string)
  
  # Get classified raster 
  x = object@timeseries$Class
  x_twdtw = object@timeseries$Distance
  
  # Reproject points to raster projection 
  y = spTransform(y, CRS(projection(x)))
  
  # Get time intervals 
  timeline = index(object)
  timeline = c(timeline[1] - diff(timeline[1:2]) - 1, timeline)
  r_intervals = data.frame(from=timeline[-length(timeline)], to=timeline[-1])
  
  # Get land use/cover classes   
  rnames = labels(object)
  rlevels = levels(object)
  
  # Compute mapped area of each class by classification interval 
  a_by_interval = lapply(1:nlayers(x), FUN = .getAreaByClass, x, rlevels, rnames)
  
  # Compute total mapped area by class 
  area_by_class = do.call("rbind", a_by_interval)
  area_by_class = colSums(area_by_class)
  
  # Get classified and predicted land cover/use classes for each control point 
  pred_classes = extract(x, y)
  pred_distance = extract(x_twdtw, y)
  samples_by_period = lapply(1:nrow(r_intervals), FUN = .getPredRefClasses, r_intervals, pred_classes, pred_distance, y, rlevels, rnames)
  samples_all = do.call("rbind", samples_by_period)
  
  # Compute error matrix 
  error_matrix_by_period = lapply(1:nrow(r_intervals), function(i) table(samples_by_period[[i]][,c("Predicted","Reference")]))
  error_matrix_summary = table(samples_all[,c("Predicted","Reference")])
  
  # Compute accuracy assessment 
  accuracy_by_period = lapply(seq_along(error_matrix_by_period), function(i) .twdtwAssess(x = error_matrix_by_period[[i]], a_by_interval[[i]], conf.int=conf.int))
  names(accuracy_by_period) = index(object)
  accuracy_summary = .twdtwAssess(error_matrix_summary, area_by_class, conf.int=conf.int)
  
  new("twdtwAssessment", accuracySummary = accuracy_summary, 
                         accuracyByPeriod = accuracy_by_period, 
                         data = samples_all)
  
}

.twdtwAssess = function(x, mapped_area, conf.int){
  
  mult = qnorm(1-(1-conf.int)/2, mean = 0, sd = 1)
  
  cnames = names(mapped_area)
  x = data.frame(cbind(x), row.names = cnames)
  cnames = names(x)
  rownames(x) = cnames
  names(mapped_area) = cnames
  
  total_map = rowSums(x)
  total_ref = colSums(x)
  total_area = sum(mapped_area)
  total_samples = sum(total_ref)
  
  # Weight 
  w = mapped_area / total_area
  
  # Error matrix 
  error_matrix = data.frame(cbind(x, Total=total_map, Area=mapped_area, w=w))
  error_matrix["Total",] = colSums(error_matrix)
  
  # Proportions 
  y = t(apply(error_matrix[!rownames(error_matrix)%in%"Total",], 1, function(x) (x[cnames] / x["Total"]) * x["w"]))
  y[total_map==0,] = 0 
  total_prop_map = rowSums(y, na.rm = TRUE)
  total_prop_ref = colSums(y, na.rm = TRUE)
  
  # Proportions matrix 
  prop_matrix = data.frame(y, Total = total_prop_map, Area = mapped_area, w = w)
  prop_matrix["Total",] = colSums(prop_matrix, na.rm = TRUE)
  
  # Accuracy 
  UA = as.numeric(diag(as.matrix(prop_matrix[cnames,cnames])) / prop_matrix[cnames,"Total"])
  UA[total_map==0] = 1 
  names(UA) = cnames
  
  PA = as.numeric(diag(as.matrix(prop_matrix[cnames,cnames])) / prop_matrix["Total",cnames])
  PA[total_ref==0] = 1 
  names(PA) = cnames
  OA = sum(diag(as.matrix(prop_matrix[cnames,cnames])), na.rm = TRUE)
  
  temp = w^2*UA*(1-UA)/(total_map-1)
  
  VO = sum(temp, na.rm = TRUE)
  SO = sqrt(VO)
  OCI = SO * mult
  
  VU = UA*(1-UA)/(total_map-1)
  SU = sqrt(VU)
  UCI = SU * mult
  
  fun1 = function(x, xt, Area){
    sum(Area*x/xt, na.rm = TRUE)
  }
  
  fun2 = function(i, x, xt, Area, PA){
    x  =  as.numeric(x[,i])
    x  =  x[-i]
    xt =  xt[-i]
    Area  =  Area[-i]
    PA = PA[i]
    PA^2*sum(Area^2*x/xt*(1-x/xt)/(xt-1), na.rm = TRUE)
  }
  
  Nj = apply(x, 2, fun1, total_map, mapped_area)
  expr1 = mapped_area^2*(1-PA)^2*UA*(1-UA)/(total_map-1)
  expr2 = sapply(1:nrow(x), fun2, x=x, xt=total_map, Area=mapped_area, PA=PA)
  VP = (1/sapply(Nj, function(x) ifelse(x==0, 1, x))^2)*(expr1+expr2)
  SP = sapply(VP, function(x) ifelse(x==0, 0, sqrt(x)))
  PCI = SP * mult
  
  # Compute adjusted area 
  estimated_area = prop_matrix["Total",cnames] * prop_matrix["Total","Area"]
  sd_error = apply(prop_matrix[cnames,cnames], 2, function(x) sqrt(sum( (prop_matrix[cnames,"w"]*x[cnames]-x[cnames]^2)/(error_matrix[cnames,"Total"]-1) )) ) 
  sd_error_estimated_area = sd_error * prop_matrix["Total","Area"]
  CI_estimated_area = sd_error_estimated_area * mult
  
  res = list(OverallAccuracy   = c(Accuracy=OA, Var=VO, sd=SO, ci=OCI),
             UsersAccuracy     = cbind(Accuracy=UA, Var=VU, sd=SU, ci=UCI),
             ProducersAccuracy = cbind(Accuracy=PA, Var=VP, sd=SP, ci=PCI),
             AreaUncertainty   = cbind(Mapped=c(prop_matrix[cnames,"Area"]), 
                                       Adjusted=c(estimated_area), 
                                       ci=c(CI_estimated_area)),
             ErrorMatrix = error_matrix,
             ProportionMatrix = prop_matrix,
             conf.int = conf.int
  )
  
  res
  
}





