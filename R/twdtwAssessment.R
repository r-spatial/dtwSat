
setGeneric("twdtwAssessment", 
           def = function(object, ...) standardGeneric("twdtwAssessment")
)

#' @inheritParams twdtwAssessment-class
#' @aliases twdtwAssessment
#' 
#' @describeIn twdtwAssessment this function performs an accuracy assessment 
#' of the classified maps. The function returns Overall Accuracy, 
#' User's Accuracy, Produce's Accuracy, and error matrix (confusion matrix) for 
#' each time interval and a summary considering all classified intervals. 
#'
#' @examples 
#' \dontrun{
#' 
#' # Create raster time series
#' 
#' evi = brick(system.file("lucc_MT/data/evi.tif", package="dtwSat"))
#' ndvi = brick(system.file("lucc_MT/data/ndvi.tif", package="dtwSat"))
#' red = brick(system.file("lucc_MT/data/red.tif", package="dtwSat"))
#' blue = brick(system.file("lucc_MT/data/blue.tif", package="dtwSat"))
#' nir = brick(system.file("lucc_MT/data/nir.tif", package="dtwSat"))
#' mir = brick(system.file("lucc_MT/data/mir.tif", package="dtwSat"))
#' doy = brick(system.file("lucc_MT/data/doy.tif", package="dtwSat"))
#' timeline = scan(system.file("lucc_MT/data/timeline", package="dtwSat"), what="date")
#' 
#' rts = twdtwRaster(evi, ndvi, red, blue, nir, mir, timeline = timeline, doy = doy)
#' field_samples = read.csv(system.file("lucc_MT/data/samples.csv", package="dtwSat"))
#' proj_str = scan(system.file("lucc_MT/data/samples_projection", 
#'                 package="dtwSat"), what = "character")
#' field_samples_ts = getTimeSeries(rts, y = field_samples, proj4string = proj_str)
#' temporal_patterns = createPatterns(field_samples_ts, freq = 8, formula = y ~ s(x))
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' 
#' # Run TWDTW analysis for raster time series 
#' 
#' r_twdtw = twdtwApply(x=rts, y=temporal_patterns, weight.fun=log_fun, format="GTiff", 
#'                      overwrite=TRUE, chunk.size=1000)
#' 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' log_fun = weight.fun=logisticWeight(-0.1,50)
#' time_interval = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), 
#'                     by="12 month")
#' r_twdtw = twdtwApply(x=rts, y=patt, weight.fun=log_fun, breaks=time_interval, 
#'           filepath="~/test_twdtw", overwrite=TRUE, format="GTiff")
#' 
#' r_lucc = twdtwClassify(r_twdtw, format="GTiff", overwrite=TRUE)
#' 
#' plotMaps(r_lucc)
#' 
#' # Map assessment 
#' 
#' 
#'   
#' }
#' @export
setMethod(f = "twdtwAssessment", 
          definition = function(object, y, id.labels, proj4string, conf.int) 
              twdtwAssessment.twdtwRaster(object, y, id.labels, proj4string, conf.int))

twdtwAssessment.twdtwRaster = function(object, y, id.labels, proj4string, conf.int){

  # Get classified raster 
  x = object@timeseries$Class
  
  # Get time intervals 
  timeline = index(object)
  timeline = c(timeline[1] - diff(timeline[1:2]) - 1, timeline)
  r_intervals = data.frame(from=timeline[-length(timeline)], to=timeline[-1])
  
  # Get land use/cover classes   
  rnames = labels(object)
  rlevels = levels(object)
  
  # Compute area of each class by classification interval 
  a_by_interval = lapply(1:nlayers(x), FUN = .area_by_class, x, rlevels, rnames)
  
  # Compute total area by class 
  area_by_class = do.call("rbind", a_by_interval)
  area_by_class = colSums(area_by_class)
  
  # Get classified and predicted land cover/use classes for each control point 
  pred_classes = extract(x, y)
  samples_by_year = lapply(1:nrow(r_intervals), FUN = .get_pre_ref_classes, r_intervals, pred_classes, y, rlevels, rnames)
  samples_all = do.call("rbind", samples_by_year)
  
  # Compute error matrix 
  error_matrix_by_year = lapply(samples_by_year, table)
  error_matrix_summary = table(samples_all)
  
  # Compute accuracy assessment 
  accuracy_by_year = lapply(seq_along(error_matrix_by_year), function(i) .twdtwAssessment(x = error_matrix_by_year[[i]], a_by_interval[[i]], conf.int))
  accuracy_summary = .twdtwAssessment(error_matrix_summary, area_by_class, conf.int=1.96)  
  
  # new("twdtwCrossValidation", partitions=partitions, accuracy=res)
  
  list(accuracy_summary, accuracy_by_year)
  
}

.get_pre_ref_classes = function(i, r_intervals, pred, y, rlevels, rnames){
  I = which((r_intervals$to[i] - as.Date(y$from) > 30) & (as.Date(y$to) - r_intervals$from[i] > 30) )
  if(length(I)<1)
    return(NULL)
  J = match(pred[I,i], rlevels)
  Predicted = factor(as.character(rnames[J]), levels = rnames, labels = rnames)
  Reference = factor(as.character(y$label[I]), levels = rnames, labels = rnames)
  data.frame(Predicted, Reference)
}

.area_by_class = function(l, r, rlevels, rnames){
  r = raster(r, layer = l)
  a = zonal(r, r, 'count')
  I = match(a[,'zone'], rlevels)
  out = rep(0, length(rnames))
  names(out) = rnames
  out[I] = a[,'count'] * prod(res(r))
  names(out) = rnames
  out  
}

.twdtwAssessment = function(x, area, conf.int){
  
  cnames = names(area)
  # cnames = paste0("aux_classname_",seq_along(cnames))
  x = data.frame(cbind(x), row.names = cnames)
  # names(x) = cnames
  cnames = names(x)
  rownames(x) = cnames
  names(area) = cnames
  
  total_map = rowSums(x)
  total_ref = colSums(x)
  total_area = sum(area)
  total_samples = sum(total_ref)
  
  # Weight 
  w = area / total_area
  
  # Error matrix 
  
  error_matrix = data.frame(cbind(x, Total=total_map, A=area, w=w))
  error_matrix["Total",] = colSums(error_matrix)
  
  # Proportions 
  y = t(apply(error_matrix[!rownames(error_matrix)%in%"Total",], 1, function(x) (x[cnames] / x["Total"]) * x["w"]))
  y[total_map==0,] = 0 
  total_prop_map = rowSums(y, na.rm = TRUE)
  total_prop_ref = colSums(y, na.rm = TRUE)
  
  # Proportions matrix 
  prop_matrix = data.frame(y, Total = total_prop_map, A = area, w = w)
  prop_matrix["Total",] = colSums(prop_matrix, na.rm = TRUE)
  
  # Accuracy 
  UA = as.numeric(diag(as.matrix(prop_matrix[cnames,cnames])) / prop_matrix[cnames,"Total"])
  UA[total_map==0] = 1 
  names(UA) = cnames
  
  PA = as.numeric(diag(as.matrix(prop_matrix[cnames,cnames])) / prop_matrix["Total",cnames])
  PA[total_ref==0] = 1 
  names(PA) = cnames
  OA = sum(diag(as.matrix(prop_matrix[cnames,cnames])), na.rm = TRUE)
  
  #a_pixel = as.numeric(prop_matrix["Total",cnames] * prop_matrix["Total","A"])
  #names(a_pixel) = cnames
  #a_ha = a_pixel*res^2/100^2
  conf.int = 1.96
  temp = w^2*UA*(1-UA)/(total_map-1)
  
  VO = sum(temp, na.rm = TRUE)
  SO = sqrt(VO)
  OCI = SO * conf.int
  
  VU = UA*(1-UA)/(total_map-1)
  SU = sqrt(VU)
  UCI = SU * conf.int
  
  fun1 = function(x, xt, A){
    sum(A*x/xt, na.rm = TRUE)
  }
  
  fun2 = function(i, x, xt, A, PA){
    x  =  as.numeric(x[,i])
    x  =  x[-i]
    xt =  xt[-i]
    A  =  A[-i]
    PA = PA[i]
    PA^2*sum(A^2*x/xt*(1-x/xt)/(xt-1), na.rm = TRUE)
  }
  
  Nj = apply(x, 2, fun1, total_map, area)
  expr1 = area^2*(1-PA)^2*UA*(1-UA)/(total_map-1)
  expr2 = sapply(1:nrow(x), fun2, x=x, xt=total_map, A=area, PA=PA)
  # VP = (1/Nj^2)*(expr1+expr2)
  VP = (1/sapply(Nj, function(x) ifelse(x==0, 1, x))^2)*(expr1+expr2)
  SP = sapply(VP, function(x) ifelse(x==0, 0, sqrt(x)))
  PCI = SP * 1.96
  
  res = list(OverallAccuracy   = c(Accuracy=OA, Var=VO, sd=SO, ci95=OCI),
             UsersAccuracy     = cbind(Accuracy=UA, Var=VU, sd=SU, ci95=UCI),
             ProducersAccuracy = cbind(Accuracy=PA, Var=VP, sd=SP, ci95=PCI),
             ErrorMatrix = error_matrix
  )
  
  res
  
}





