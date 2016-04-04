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
#' patt = twdtwTimeSeries(patterns.list)
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


#' @title Classification assessment 
#' @name Assessment
#' 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This functions create data partitions and compute assessment metrics. 
#' 
#' @param object an object of class \code{\link[dtwSat]{twdtwTimeSeries}} or 
#' \code{\link[dtwSat]{twdtwMatches}}.
#' 
#' @param times Number of partitions to create.
#' 
#' @param p the percentage of data that goes to training. 
#' See \code{\link[caret]{createDataPartition}} for details.
#' 
#' @param ... Other arguments to be passed to \code{\link[dtwSat]{createPatterns}}.
#' 
#' @param matrix logical. If TRUE retrieves the confusion matrix. 
#' FALSE retrieves User's Accuracy (UA) and Producer's Accuracy (PA). 
#' Dafault is FALSE. 
#' 
#' @details 
#' \describe{
#'  \item{\code{splitDataset}:}{This function splits the a set of time 
#'        series into training and validation. The function uses stratified 
#'        sampling and a simple random sampling for each stratum. Each data partition 
#'        returned by this function has the temporal patterns and a set of time series for 
#'        validation.}
#'  \item{\code{twdtwAssess}:}{The function \code{splitDataset} performs the assessment of 
#'        the classification based on the labels of the classified time series 
#'        (Reference) and the labels of the classification (Predicted). This function
#'        returns a data.frame with User's and Produce's Accuracy or a list for confusion 
#'        matrices.}
#' }
#'
#' @seealso 
#' \code{\link[dtwSat]{twdtwMatches-class}},
#' \code{\link[dtwSat]{twdtwApply}}, and 
#' \code{\link[dtwSat]{twdtwClassify}}.
#'
#' @examples
#' \dontrun{
#' load(system.file("lucc_MT/field_samples_ts.RData", package="dtwSat"))
#' set.seed(1)
#' partitions = splitDataset(field_samples_ts, p=0.1, times=5, 
#'                           freq = 8, formula = y ~ s(x, bs="cc"))
#' log_fun = logisticWeight(alpha=-0.1, beta=50) 
#' twdtw_res = lapply(partitions, function(x){
#'    res = twdtwApply(x = x$ts, y = x$patterns, weight.fun = log_fun, n=1)
#'    twdtwClassify(x = res)
#' })
#' assessment = twdtwAssess(twdtw_res)
#' head(assessment, 5)
#' }
NULL

setGeneric("splitDataset", function(object, times, p, ...) standardGeneric("splitDataset"))

#' @rdname Assessment
#' @aliases splitDataset
#' @export
setMethod("splitDataset", "twdtwTimeSeries",
          function(object, times=1, p=0.1, ...) splitDataset.twdtwTimeSeries(object, times=times, p=p, ...))
            
splitDataset.twdtwTimeSeries = function(object, times, p, ...){
  
  partitions = createDataPartition(y = labels(object), times, p, list = TRUE)
  
  res = lapply(partitions, function(I){
      training_ts = subset(object, I)
      validation_ts = subset(object, -I)
      patt = createPatterns(training_ts, ...)
      list(patterns=patt, ts=validation_ts)
  })
  
  res
}


setGeneric("twdtwAssess", function(object, matrix=FALSE) standardGeneric("twdtwAssess"))

#' @rdname Assessment
#' @aliases twdtwAssess
#' @export
setMethod("twdtwAssess", "list",
          function(object, matrix) twdtwAssess.twdtwTimeSeries(object, matrix=matrix))

twdtwAssess.twdtwTimeSeries = function(object, matrix){

 res = lapply(object, function(x){
        ref = labels(x)$timeseries
        levels = sort(as.character(unique(ref)))
        labels = levels 
        pred = factor(do.call("rbind", x[])$label, levels, labels)
        ref = factor(ref, levels, labels)
        table(Reference=ref, Predicted=pred)
 })
    
 if(!matrix){
    res = do.call("rbind", lapply(seq_along(res), function(i){
        x = res[[i]]
        Users = diag(x) / rowSums(x)
        Producers = diag(x) / colSums(x)
        data.frame(resample=i,label=names(Users), UA = Users, PA = Producers, row.names=NULL)
    }))
  }
  
  res
  
}


#' @title Create time sequence
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a sequence of dates for 
#' each year. The sequences start on January 1st of each year.
#' 
#' @param year A vector with the years. Default 
#' is form 2000 to the system time year \code{format(Sys.time(), ''\%Y'')}
#' @param frequency An integer with the frequency in days. Default is 16 days
#' @docType methods
#' 
#' @return vector of class \code{\link[base]{Dates}} 
#' 
#' @seealso \link[dtwSat]{twdtwTimeSeriesFromVector} 
#' 
#' @examples
#' dates = createTimeSequence()
#' dates
#' 
#' @export
createTimeSequence = function(year=2000:format(Sys.time(), "%Y"), frequency=16){
  res = unlist(lapply(year, function(y){
    days = seq(from = as.Date(paste0(y,"-01-01")), to = as.Date(paste0(y,"-12-31")), by = frequency)
  }))
  res = as.Date(res, origin="1970-01-01")
  res
}


#' @title Create time series from vectors 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function creates a twdtwTimeSeries object using data stored in 
#' vectors.
#' 
#' @param x A vector of integers with columns indices. 
#' @param y A vector of integers with rows indices. 
#' @param timeline A vector of dates with the date in the format ''YYYY-MM-DD''. 
#' @param ... numeric vectors with the data. This are the attributes of the time series.  
#' 
#' @details All arguments passed to this function must have the same length. 
#' 
#' @return An object of class \code{\link[dtwSat]{twdtwTimeSeries}} 
#' 
#' @seealso \link[dtwSat]{createTimeSequence} 
#' 
#' @examples
#' n = seq(100)
#' col_id = rep(c(1,2), each=length(n))
#' row_id = rep(1, length(col_id))
#' from = as.Date("2000-01-01")
#' dates = rep(seq(from, from+length(n)-1, by=1), 2)
#' 
#' t = seq(0,4*pi,,length(col_id))
#' attr1 = 3*sin(2*t)+runif(seq(length(col_id)))*2
#' attr2 = 3*sin(2*t)+rnorm(seq(length(col_id)))*2
#' 
#' ts = twdtwTimeSeriesFromVector(x=col_id, y=row_id, timeline=dates, attr1, attr2)
#' plot(ts, type="timeseries")
#' 
#' @export
twdtwTimeSeriesFromVector = function(x, y, timeline, ...){
  arg_names = names(list(...))
  not_named = setdiff(as.character(match.call(expand.dots=TRUE)), as.character(match.call(expand.dots=FALSE)))
  if(is.null(arg_names)){ 
    arg_names = not_named
  } else {
    arg_names[arg_names==""] = not_named[arg_names==""]
  }
  data_list = list(...)
  names(data_list) = c(arg_names)
  
  I = unique(y)
  J = unique(x)
  indexArray = list()
  k = 1
  for(i in I)
    for(j in J)
    {
      indexArray[[k]] = c(i=i, j=j)
      k = k + 1
    }
  
  ts_zoo = lapply(indexArray, function(p){
    idx = which(x==p["j"] & y==p["i"])
    idx = idx[!duplicated(timeline[idx])] # Remove duplicated dates 
    if( length(idx) < 2 )
      return(NULL)
    
    datasets = lapply(data_list, function(data){
      data[idx]
    })
    zoo(data.frame(datasets), order.by = timeline[idx])
  })
  
  labels = sapply(indexArray, paste, collapse="_")
  
  twdtwTimeSeries(ts_zoo, labels = labels)
  
}



#' @title Get MODIS col row from longitude and latitude 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the MODIS col row  
#' from pair of coordinates.
#' 
#' @param x An object of class SpatialPoints or data.frame with 
#' longitude and latitude.
#' @param pixelsize A number. The MODIS resolution.
#' @param proj4string A CRS object. See \link[sp]{CRS} for details. 
#' 
#' @seealso \link[dtwSat]{modisColRowToLongLat} 
#' 
#' @examples 
#' # Example for MOD13Q1 spatial resolution 231.6564 m 
#' coord = data.frame(-55.93715, -12.05729)
#' prj_str = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#' modisColRowFromLongLat(x=coord, pixelsize=231.6564, proj4string=prj_str)
#' 
#' @export
modisColRowFromLongLat = function(x, pixelsize, proj4string=CRS(as.character(NA))){
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  tilesize = 1111950.51966666709631681442
  nhtiles = 36
  nvtiles = 18
  ncells = trunc(tilesize/pixelsize)
  xmin = - tilesize*nhtiles/2
  ymax = tilesize*nvtiles/2
  
  if(is(x, "data.frame")){
    if(!is(proj4string, "CRS")) proj4string = try(CRS(proj4string))
    I = match(c("longitude","latitude"), names(x))
    if(any(is.na(I))|length(I)<2) I = c(1,2)
    sppoint = SpatialPoints(x[,I], proj4string = proj4string)
  }

  longitude = coordinates(sppoint)[,1]
  latitude = coordinates(sppoint)[,2]
  sppoint = spTransform(sppoint, projCRSSinu)
  x = coordinates(sppoint)[,1]
  y = coordinates(sppoint)[,2]
  
  dx = x - xmin
  dy = ymax - y 
  
  col = trunc(dx/pixelsize)
  row = trunc(dy/pixelsize)
  
  h = abs(trunc(dx / tilesize))
  v = abs(trunc(dy / tilesize))
  
  j = col - ncells * h
  i = row - ncells * v
  
  res = data.frame(longitude, latitude, col, row, h, v, j, i)
  rownames(res) = NULL
  res
}


#' @title Get longitude and latitude from MODIS col row from 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function retrieves the longitude and latitude from 
#' MODIS col row.
#' 
#' @param col An integer. A global index of the MODIS column
#' or a local index of the column within a specific h:v tile.
#' @param row An integer. A global index of the MODIS row
#' or a local index of the row within a specific h:v tile.
#' @param h An integer. Horizontal MODIS tile.
#' @param v An integer. Vertical MODIS tile.
#' @param pixelsize A number. The MODIS resolution.
#' @param proj4string A CRS object. See \link[sp]{CRS} for details. 
#' 
#' 
#' @seealso \link[dtwSat]{modisColRowFromLongLat} 
#' 
#' @examples 
#' # Example for MOD13Q1 spatial resolution 231.6564 m  
#' modisColRowToLongLat(col=60124, row=48978, pixelsize=231.6564)
#' 
#' @export
modisColRowToLongLat = function(col, row, h=NULL, v=NULL, pixelsize, proj4string=NULL){
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  
  if(is(proj4string, "character")) proj4string = try(CRS(proj4string))
  if(!is(proj4string,"CRS")) proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  tilesize = 1111950.51966666709631681442 # meters
  ncells = trunc(tilesize/pixelsize)
  nhtiles = 36
  nvtiles = 18
  
  nhcells = ncells*nhtiles/2
  nvcells = ncells*nvtiles/2
  j = col - nhcells
  i = nvcells - row
  if( !is.null(h) )
    j = ncells*(h - nhtiles/2) + col
  if( !is.null(v) )
    i = ncells*(nvtiles/2 - v) - row 
  
  longitude = pixelsize*j + pixelsize/2
  latitude = pixelsize*i - pixelsize/2
  sppoint = SpatialPoints(coords = data.frame(longitude, latitude), projCRSSinu)
  sppoint = spTransform(sppoint, proj4string)
  res = data.frame(coordinates(sppoint))
  res
}


