
#' @title Get dates from year and day of the year
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the date corresponding to the given 
#' year and day of the year.
#' 
#' @param year A vector with the years.
#' @param doy A vector with the day of the year. 
#' It must have the same length as \code{year}.
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{Dates}} object.
#' 
#' @seealso \link[dtwSat]{shiftDates} 
#' 
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples
#' year = c(2000, 2001)
#' doy = c(366, 365)
#' dates = getDatesFromDOY(year, doy)
#' dates
#'
#' @export
getDatesFromDOY = function(year, doy){
  res = as.Date(paste(as.numeric(year), as.numeric(doy)), format="%Y %j", origin="1970-01-01")
  I = which(diff(res)<0)+1
  res[I] = as.Date(paste0(as.numeric(format(res[I],"%Y"))+1, format(res[I], "-%m-%d")))
  res
}



#' @title Shift dates 
#' @name shiftDates
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function shifts the dates of the time series to a 
#' given base year. 
#' 
#' @param object \code{\link[dtwSat]{twdtwTimeSeries}} objects, 
#' \code{\link[zoo]{zoo}} objects or a list of \code{\link[zoo]{zoo}} objects.
#' 
#' @param year the base year to shift the time series to. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}
#'
#' @return An object of the same class as the input \code{object}. 
#'
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @export
setGeneric("shiftDates", function(object, year=NULL) standardGeneric("shiftDates"))

#' @rdname shiftDates
#' @aliases shiftDates-twdtwTimeSeries
#' 
#' @references
#'   \insertRef{Maus:2019}{dtwSat}
#'   
#'   \insertRef{Maus:2016}{dtwSat}
#'   
#' @examples
#' patt = twdtwTimeSeries(MOD13Q1.patterns.list)
#' npatt = shiftDates(patt, year=2005)
#' index(patt)
#' index(npatt)
#' 
#' @export
setMethod("shiftDates", "twdtwTimeSeries",
          function(object, year) 
            do.call("twdtwTimeSeries", lapply(as.list(object), FUN=shiftDates.twdtwTimeSeries, year=year)))

#' @rdname shiftDates
#' @aliases shiftDates-list
#' @export
setMethod("shiftDates", "list",
          function(object, year) 
            shiftDates(twdtwTimeSeries(object), year=year)[])

setOldClass("zoo")
#' @rdname shiftDates
#' @aliases shiftDates-zoo
#' @export
setMethod("shiftDates", "zoo",
          function(object, year) 
            shiftDates(twdtwTimeSeries(object), year=year)[[1]])

            
shiftDates.twdtwTimeSeries = function(x, year){
  labels = as.character(labels(x))
  x = x[[1]]
  dates = index(x)
  last_date = tail(dates, 1)
  shift_days = as.numeric(last_date - as.Date(paste0(year,format(last_date, "-%m-%d"))))
  d = as.numeric(dates) - shift_days
  new("twdtwTimeSeries", timeseries=zoo(data.frame(x), as.Date(d)), labels=labels)
}


.adjustFactores = function(ref, pred, levels=NULL, labels=NULL){
  ref  = as.character(ref)
  pred = as.character(pred)
  if(is.null(levels))
    levels = sort(unique(ref))
  if(is.null(labels))
    labels = levels
  ref  = factor(ref,  levels, labels)
  pred = factor(pred, levels, labels)
  data = data.frame(Predicted=pred, Reference=ref)
}

.adjustLabelID = function(y, labels, id.labels){
  if(!"label"%in%names(y)) y$label = paste0("ts",row.names(y))
  if(!is.null(id.labels)) y$label = as.character(y[[id.labels]])
  if(!is.null(id.labels) & !is.null(labels)){
    I = which(!is.na(match(as.character(y$label), as.character(labels))))
    if(length(I)<1) 
      stop("There are no matches between id.labels and labels")
  } else if(!is.null(labels)) { 
    y$label = as.character(labels)
  }
  y
}

.toSpatialPointsDataFrame = function(y, object, proj4string){
  if(is(y, "data.frame")){
    if(is.null(proj4string)){
      warning("Missing projection. Setting the same projection as the raster time series.", call. = FALSE)
      proj4string = CRS(projection(object@timeseries[[1]]))
    }
    if(!is(proj4string, "CRS")) proj4string = try(CRS(proj4string))
    y = SpatialPointsDataFrame(y[,c("longitude","latitude")], y, proj4string = proj4string)
  }
  if(!(is(y, "SpatialPoints") | is(y, "SpatialPointsDataFrame")))
    stop("y is not SpatialPoints or SpatialPointsDataFrame")
  row.names(y) = 1:nrow(y)
  y
}


.getPredRefClasses = function(i, r_intervals, pred_classes, pred_distance, y, rlevels, rnames){
  i_leng = as.numeric(r_intervals$to[i] - r_intervals$from[i])
  from   = as.Date(y$from)
  to     = as.Date(y$to)
  # Select overlapping alignments 
  J      = which(from <= r_intervals$to[i] & to >= r_intervals$from[i])
  # Adjust overlapping 
  from   = sapply(from[J], function(x) ifelse(x < r_intervals$from[i], r_intervals$from[i], x))
  to     = sapply(to[J], function(x) ifelse(x > r_intervals$to[i], r_intervals$to[i], x))
  # Compute overlapping proportion 
  if(length(to)<1)
    return(NULL)
  i_over = to - from 
  # print(i_leng)
  # print(i_over)
  prop_over = abs(i_over / i_leng)
  # Select alignments 
  I = which(prop_over > .5)
  # I = which((r_intervals$to[i] - as.Date(y$from) > 30) & (as.Date(y$to) - r_intervals$from[i] > 30) )
  if(length(J[I])<1)
    return(NULL)
  K = match(pred_classes[J[I],i], rlevels)
  Predicted = factor(as.character(rnames[K]), levels = rnames, labels = rnames)
  Reference = factor(as.character(y$label[J[I]]), levels = rnames, labels = rnames)
  Distance = pred_distance[J[I],i]
  data.frame(Sample.id = row.names(y)[J[I]], coordinates(y[J[I],]), Period=i, from=r_intervals$from[i], to=r_intervals$to[i], Predicted, Reference, Distance)
}

.getAreaByClass = function(l, r, rlevels, rnames){
  r = raster(r, layer = l)
  if(isLonLat(r)){
    warning("Computing the approximate surface area in km2 of cells in an unprojected (longitude/latitude) Raster object. See ?raster::area", call. = TRUE)
    # r = projectRaster(from = r, crs = proj_str, method = 'ngb')
    ra = area(r)
    I = lapply(rlevels, function(i) r[]==i )
    out = sapply(I, function(i) sum(ra[i], na.rm = TRUE) )
    names(out) = rnames
  } else {
    npx = zonal(r, r, 'count')
    I = match(npx[,'zone'], rlevels)
    out = rep(0, length(rnames))
    names(out) = rnames
    out[I] = npx[,'count'] * prod(res(r))
    names(out) = rnames
  }
  out
}


