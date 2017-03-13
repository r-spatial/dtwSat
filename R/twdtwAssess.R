
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
#' is an error matrix (confusion matrix). If \code{object} is \code{\link[dtwSat]{twdtwMatches}} 
#' area can be either a vector with the area of each classified object, or a single number 
#' if the objects are single pixels. 

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
#' @param rm.nosample if sum of columns and sum of rows of the error matrix are zero 
#' then remove class. Default is TRUE. 
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
NULL

#' @aliases twdtwAssess-twdtwRaster
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
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
#' # Get time series form raster
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' validation_ts = getTimeSeries(rts, y = validation_samples, proj4string = proj_str)
#' 
#' # Create temporal patterns 
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
#'                            proj4string = proj_str, conf.int = .95, rm.nosample=TRUE) 
#' twdtw_assess
#' 
#' # Plot assessment 
#' plot(twdtw_assess, type="accuracy")
#' plot(twdtw_assess, type="area")
#' plot(twdtw_assess, type="map", samples = "all")
#' plot(twdtw_assess, type="map", samples = "incorrect")
#' plot(twdtw_assess, type="map", samples = "correct")
#' 
#' # Create latex tables 
#' twdtwXtable(twdtw_assess, table.type="matrix")
#' twdtwXtable(twdtw_assess, table.type="accuracy")
#' twdtwXtable(twdtw_assess, table.type="area")
#' 
#' }
#' @export
setMethod(f = "twdtwAssess", signature = "twdtwRaster",
          definition = function(object, y, labels=NULL, id.labels=NULL, proj4string=NULL, conf.int=.95, rm.nosample=TRUE) 
            twdtwAssess.twdtwRaster(object, y, labels, id.labels, proj4string, conf.int, rm.nosample))

#' @aliases twdtwAssess-data.frame
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' @export
setMethod(f = "twdtwAssess", signature = "data.frame",
          definition = function(object, area, conf.int=.95, rm.nosample=TRUE) 
            twdtwAssess.table(object, area, conf.int, rm.nosample))

#' @aliases twdtwAssess-table
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' @export
setMethod(f = "twdtwAssess", signature = "table",
          definition = function(object, area, conf.int=.95, rm.nosample=TRUE) 
            twdtwAssess(as.data.frame.matrix(object), area, conf.int, rm.nosample))

#' @aliases twdtwAssess-matrix
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
#' 
#' @examples 
#' 
#' # Total mapped area by class. Data from [1]
#' area = c(A = 22353, B = 1122543, C = 610228) 
#' 
#' # Error matrix, columns (Reference) rows (Map)
#' x = 
#'     rbind(
#'          c( 97,	 0,   3),
#'          c(  3, 279,  18),
#'          c(  2,   1,  97)
#'    )
#'
#' table_assess = twdtwAssess(x, area, conf.int = .95)
#' 
#' table_assess
#' 
#' plot(table_assess, type="area", perc=FALSE)
#' 
#' plot(table_assess, type="accuracy")
#' 
#' @export
setMethod(f = "twdtwAssess", signature = "matrix",
          definition = function(object, area, conf.int=.95, rm.nosample=TRUE) 
            twdtwAssess(as.data.frame.matrix(object), area, conf.int, rm.nosample))

#' @aliases twdtwAssess-twdtwMatches
#' @inheritParams twdtwAssess
#' @rdname twdtwAssess 
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
#' # Get time series form raster
#' training_ts = getTimeSeries(rts, y = training_samples, proj4string = proj_str)
#' validation_ts = getTimeSeries(rts, y = validation_samples, proj4string = proj_str)
#' 
#' # Create temporal patterns 
#' temporal_patterns = createPatterns(training_ts, freq = 8, formula = y ~ s(x))
#' 
#' # Run TWDTW analysis for raster time series 
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' r_twdtw = twdtwApply(x = validation_ts, 
#'                      y = temporal_patterns, weight.fun = log_fun)
#' 
#' # Accuracy assessment 
#' twdtw_assess = twdtwAssess(r_twdtw, area = 53664.67, conf.int=.95)
#' twdtw_assess
#' 
#' plot(twdtw_assess, type="accuracy")
#' plot(twdtw_assess, type="area")
#' 
#' twdtwXtable(twdtw_assess, table.type="matrix")
#' twdtwXtable(twdtw_assess, table.type="accuracy")
#' twdtwXtable(twdtw_assess, table.type="area")
#' 
#' }
#' @export
setMethod(f = "twdtwAssess", signature = "twdtwMatches",
          definition = function(object, area, conf.int, rm.nosample) 
            twdtwAssess.twdtwTimeSeries(object, area, conf.int, rm.nosample))

twdtwAssess.twdtwTimeSeries = function(object, area, conf.int, rm.nosample){
  
  df = do.call("rbind", lapply(object[], function(xx) xx[which.min(xx$distance),]) )
  
  ref = labels(object)$timeseries
  
  pred = as.character(df$label)
  
  data = data.frame(.adjustFactores(ref, pred, levels=NULL, labels=NULL), df[,!names(df)%in%"labels"])

  error_matrix = table(Predicted=data$Predicted, Reference=data$Reference)
  
  if(length(area)==1)
    a = rep(area, length(object@timeseries))

  a = aggregate(x = a, by = list(pred), FUN = sum)
  
  area = a$x 
  
  names(area) = a$Group.1
  
  res = .twdtwAssess(error_matrix, area, conf.int, rm.nosample)

  new("twdtwAssessment", accuracySummary=res)
}

twdtwAssess.table = function(object, area, conf.int, rm.nosample){
  
  if(ncol(object)!=nrow(object))
    stop("object has have the same number of rows and columns")
  
  if(nrow(object)!=length(area))
    stop("area must have length equal to the number of rows in object")
  
  accuracy = .twdtwAssess(object, area, conf.int, rm.nosample)
  
  new("twdtwAssessment", accuracySummary=accuracy, accuracyByPeriod=accuracy)
}
  
twdtwAssess.twdtwRaster = function(object, y, labels, id.labels, proj4string, conf.int, rm.nosample){
  
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
  y = spTransform(y, CRS(projection(object)))
  
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
  accuracy_by_period = lapply(seq_along(error_matrix_by_period), function(i) 
    .twdtwAssess(x = error_matrix_by_period[[i]], a_by_interval[[i]], conf.int=conf.int, rm.nosample=FALSE))
  names(accuracy_by_period) = index(object)
  accuracy_summary = .twdtwAssess(error_matrix_summary, area_by_class, conf.int=conf.int, rm.nosample)
  
  sp.data = SpatialPointsDataFrame(coords = samples_all[,c("longitude", "latitude")], 
                                   data = samples_all[,!names(samples_all)%in%c("longitude", "latitude")],
                                   proj4string = CRS(projection(object)))
  
  new("twdtwAssessment", accuracySummary = accuracy_summary, 
                         accuracyByPeriod = accuracy_by_period, 
                         data = sp.data,
                         map = object)
  
}

.twdtwAssess = function(x, mapped_area, conf.int, rm.nosample){
  
  mult = qnorm(1-(1-conf.int)/2, mean = 0, sd = 1)
  
  cnames = names(mapped_area)
  rownames(x) = cnames
  names(x) = cnames
  
  total_map = rowSums(x)
  total_ref = colSums(x)
  total_area = sum(mapped_area)
  total_samples = sum(total_ref)
  
  if(rm.nosample){
    I = total_ref>0 | total_map>0
    x = x[I,I]
    cnames = cnames[I]
    total_map = total_map[I]
    total_ref = total_ref[I]
    mapped_area = mapped_area[I]
    total_area = sum(mapped_area)
    total_samples = sum(total_ref)
  }
  
  # Weight 
  w = mapped_area / total_area
  
  # Error matrix 
  error_matrix = cbind(x, Total=total_map, Area=mapped_area, w=w)
  error_matrix = rbind(error_matrix, Total = colSums(error_matrix))
  
  # Proportions 
  y = t(apply(error_matrix[!rownames(error_matrix)%in%"Total",], 1, function(x) (x[cnames] / x["Total"]) * x["w"]))
  y[total_map==0,] = 0 
  total_prop_map = rowSums(y, na.rm = TRUE)
  total_prop_ref = colSums(y, na.rm = TRUE)
  
  # Proportions matrix 
  prop_matrix = cbind(y, Total = total_prop_map, Area = mapped_area, w = w)
  prop_matrix = rbind(prop_matrix, Total = colSums(prop_matrix, na.rm = TRUE))
  
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





