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
#   R Package dtwSat - 2016-01-16                             #
#                                                             #
###############################################################


#' @title Get dates from year and day of the year
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the date corresponding to the ginven 
#' year and day of the year.
#' 
#' @param year An vector with the years.
#' @param doy An vector with the day of the year. 
#' It must have the same lenght as \code{year}.
#' 
#' @docType methods
#' 
#' @return A \code{\link[base]{Dates}} object.
#' 
#' @seealso \link[dtwSat]{shiftDates} 
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
#' @param year the base year to shit the time series. 
#' 
#' @seealso 
#' \code{\link[dtwSat]{twdtwTimeSeries-class}}
#'
#' @return An object of the same class as the input \code{object}. 
#'
#' @export
setGeneric("shiftDates", function(object, year=NULL) standardGeneric("shiftDates"))

#' @rdname shiftDates
#' @aliases shiftDates-twdtwTimeSeries
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
      stop("there is no matches between id.labels and labels")
  } else if(!is.null(labels)) { 
    y$label = as.character(labels)
  }
  y
}

.toSpatialPointsDataFrame = function(y, object, proj4string){
  if(is(y, "data.frame")){
    if(is.null(proj4string)){
      warning("Missing projection. Setting the same projection as the raster time series.", call. = FALSE)
      proj4string = CRS(projection(object))
    }
    if(!is(proj4string, "CRS")) proj4string = try(CRS(proj4string))
    y = SpatialPointsDataFrame(y[,c("longitude","latitude")], y, proj4string = proj4string)
  }
  if(!(is(y, "SpatialPoints") | is(y, "SpatialPointsDataFrame")))
    stop("y is not SpatialPoints or SpatialPointsDataFrame")
  row.names(y) = 1:nrow(y)
  y
}
