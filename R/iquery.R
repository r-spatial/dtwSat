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
#   R Package dtwSat - 2016-04-08                             #
#                                                             #
###############################################################


#' @title TWDTW analysis on SciDB 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function calls the method twdtwApply and retrieves the results 
#' as a list.
#' 
#' @inheritParams twdtwApply 
#' 
#' @details The output of this function is a list with the attributes used by r_exec to save 
#' a SciDB array. The labels of the input time series must have the format "<row>_<col>", where 
#' "row" and "col" are the indices of the SciDB array where the time sereis come from. 
#' 
#' @return An object of class \code{\link[base]{list}} whose attributes are:
#' "time_id", "class_id", "distance", "ts_id", "row_id", "col_id", "cdoy", and "year". 
#' 
#' @seealso \link[dtwSat]{twdtwApply} 
#' 
#' @examples
#' log_fun = logisticWeight(-0.1, 100)
#' ts = twdtwTimeSeries(example_ts.list, labels = paste(c(1,2), c(1,1), sep="_") )
#' patt = twdtwTimeSeries(patterns.list)
#' time_intervals = seq(from=as.Date("2007-09-01"), to=as.Date("2013-09-01"), by="6 month")
#' mat1 = twdtwApplySciDB(x=ts, y=patt, weight.fun=log_fun, breaks=time_intervals)
#' mat1
#' 
#' @export
twdtwApplySciDB = function(x, y, resample=TRUE, length=NULL, weight.fun=NULL, 
                           dist.method="Euclidean", step.matrix = symmetric1, n=NULL, 
                           span=NULL, min.length=0.5, theta = 0.5, 
                           breaks=NULL, from=NULL, to=NULL, by=NULL, overlap=0.5, ...){
                  if(is.null(breaks))
                    if( !is.null(from) &  !is.null(to) ){
                      breaks = seq(as.Date(from), as.Date(to), by=by)
                    } else {
                      patt_range = lapply(index(y), range)
                      patt_diff = trunc(sapply(patt_range, diff)/30)+1
                      min_range = which.min(patt_diff)
                      by = patt_diff[[min_range]]
                      cycles = c(18,12,6,4,3,2)
                      by = cycles[which.min(abs(by-cycles))]
                      from = patt_range[[min_range]][1]
                      to = from 
                      month(to) = month(to) + by
                      dates = as.Date(unlist(index(x)))
                      year(from) = year(min(dates))
                      year(to) = year(max(dates))
                      breaks = seq(from, to, paste(by,"month"))
                    }
                  breaks = as.Date(breaks)
                  .twdtwApply.SciDB(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, breaks, overlap, ...)
}

.twdtwApply.SciDB = function(x, y, weight.fun, dist.method, step.matrix, n, span, min.length, theta, 
                                  breaks, overlap, 
                                  mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, 
                                  mc.cores = getOption("mc.cores", 1L), mc.cleanup = TRUE, ...){
  get_aligs = function(x){
    # twdtwApply(x, y=y, weight.fun=weight.fun)
    twdtwApply(x, y=y, weight.fun=weight.fun, dist.method=dist.method, step.matrix=step.matrix, n=n, span=span, min.length=min.length, theta=theta, keep=FALSE)
  }

  # Set raster levels and labels 
  levels = levels(y)
  names(levels) = levels
  m = length(levels)
  n = length(breaks)-1
  
  # Apply TWDTW analysis  
  #twdtw_results = lapply(as.list(ts), FUN=get_aligs)
  twdtw_results = mclapply(as.list(ts), FUN=get_aligs, mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores, mc.cleanup = mc.cleanup)
  
  # Get best mathces for each point, period, and pattern 
  #res = lapply(twdtw_results, FUN=.lowestDistances, m=m, n=n, levels=levels, breaks=breaks, overlap=overlap, fill=9999)  
  res = mclapply(twdtw_results, FUN=.lowestDistances, m=m, n=n, levels=levels, breaks=breaks, overlap=overlap, fill=9999, mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores, mc.cleanup = mc.cleanup)
  # Reshape list to array 
  res = melt(res)
  colRow = lapply(as.character(labels(x)), function(a) as.numeric(unlist(strsplit(a, split = "_"))))
  res = cbind(res, do.call("rbind", colRow[res$L1]))
  names(res) = c("time_id", "class_id", "distance", "ts_id", "row_id", "col_id")
  res$cdoy = as.numeric(format(breaks[-1], "%j"))[res$time_id]
  res$year = as.numeric(format(breaks[-1], "%Y"))[res$time_id]
  lapply(c(res[c("col_id", "row_id", "time_id", "class_id", "ts_id", "year", "cdoy", "distance")]), as.double)
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



#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an aql query to create 
#' the output array for the TWDTW analysis. 
#' 
#' @param array A character for the array name. 
#' @param index.start the starting indices for each dimensions c(<col>, <row>, <time>, <class>).
#' @param index.end the ending indices for each dimensions c(<col>, <row>, <time>, <class>). 
#' @param chunk the chunk size for each dimensions c(<col>, <row>, <time>, <class>).
#' @param overlap the overlap for each dimensions c(<col>, <row>, <time>, <class>). 
#' 
#' @seealso \link[dtwSat]{modisColRowFromLongLat} 
#' 
#' @examples 
#' 
#' 
#' @export
iqueryCreateTWDTWArray = function(array, index.start, index.end, chunk, overlap){
  
  col_ids = paste(index.start[1], index.end[1], sep=":")
  row_ids = paste(index.start[2], index.end[2], sep=":")
  time_ids = paste(index.start[3], index.end[3], sep=":")
  class_ids = paste(index.start[4], index.end[4], sep=":")
  
  schema = paste("<ts_id:double,year:double,cdoy:double,distance:double>",
                 " [col_id=",col_ids,",",chunk[1],",",overlap[1],
                 ",row_id=",row_ids,",",chunk[2],",",overlap[2],
                 ",time_id=",time_ids,",",chunk[3],",",overlap[3],
                 ",class_id=",class_ids,",",chunk[4],",",overlap[4],"]", sep="")

  res = paste("CREATE ARRAY ",array, schema, sep="")
  res
}

#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an aql query to create 
#' the output array for the TWDTW analysis. 
#' 
#' @param array A character for the array name. 
#' @param index.start the starting indices for each dimensions c(<col>, <row>, <time>, <class>).
#' @param index.end the ending indices for each dimensions c(<col>, <row>, <time>, <class>). 
#' @param chunk the chunk size for each dimensions c(<col>, <row>, <time>, <class>).
#' @param overlap the overlap for each dimensions c(<col>, <row>, <time>, <class>). 
#' 
#' @seealso \link[dtwSat]{modisColRowFromLongLat} 
#' 
#' @examples 
#' 
#' str_iquery = iqueryCreateArray(array = "TEST_ARRAY",
#'                                min.col = 59354, max.col = 60132, 
#'                                min.row = 48591, max.row = 49096)
#' 
#' str_iquery
#' 
#' @export
iqueryCreateRankArray = function(array, index.start, index.end, chunk, overlap){
  
  col_ids = paste(index.start[1], index.end[1], sep=":")
  row_ids = paste(index.start[2], index.end[2], sep=":")
  time_ids = paste(index.start[3], index.end[3], sep=":")
  class_ids = paste(index.start[4], index.end[4], sep=":")
  rank_ids = paste(index.start[4], index.end[4], sep=":")
  
  schema = paste("<distance:double>",
                 " [col_id=",col_ids,",",chunk[1],",",overlap[1],
                 ",row_id=",row_ids,",",chunk[2],",",overlap[2],
                 ",time_id=",time_ids,",",chunk[3],",",overlap[3],
                 ",class_id=",class_ids,",",chunk[4],",",overlap[4],
                 ",rank_id=",rank_ids,",",chunk[4],",",overlap[4],"]", sep="")
  
  res = paste("CREATE ARRAY ",array, schema, sep="")
  res
}


#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an afl query to redimension 
#' a subset of the input array for the TWDTW analysis. 
#' 
#' @param array.in A character for the input array name.
#' @param array.out A character for the output array name. If not declared the 
#' iquery will not same the final results.
#' @param attributes A character vector for the array attributes. 
#' @param expr An R expression for processing. The expression must retrieve a 
#' list of 8 elements of type double. 
#' @param grid An object of class \link[sp]{SpatialPolygons}. Each polygon is 
#' one processing block. 
#' @param pixelsize A number. The MODIS product resolution.
#' @param blocksize processing block size c(<col>, <row>). Used instead 
#' of \code{grid} and \code{pixelsize}. 
#' @param min minimum indices for subset c(<col>, <row>, <time>).
#' @param max maximum indices for subset c(<col>, <row>, <time>).
#' @param chunk scidb chunk size for input c(<col>, <row>, <time>).
#' @param overlap scidb overlap for for input c(<col>, <row>, <time>). 
#' 
#' @seealso \link[dtwSat]{twdtwApplyiquery} 
#' 
#' @examples 
#' 
#' @export
iqueryTWDTWBlock = function(array.in, array.out, attributes, expr, grid = NULL, pixelsize = NULL,
                            blocksize, min, max, chunk, overlap=c(0,0,0)){
  
  scidbarray.in = scidb(array.in)
  scidbarray.out = scidb(array.out)
  index.start = as.numeric(scidb_coordinate_start(scidbarray.in))
  index.end = as.numeric(scidb_coordinate_end(scidbarray.in))
  outputschema = schema(scidbarray.out)
  
  if(missing(min)) min = index.start
  if(missing(max)) max = index.end
  if(length(min)<3) min[3] = index.start[3]
  if(length(max)<3) max[3] = index.end[3]
      
  if(is.null(grid)){
    blocks = .blockSizeFromColRow(min, max, blocksize)
  } else {
    blocks = .blockSizeFromPoly(min, max, grid, pixelsize)
  }

  res = lapply(blocks, function(b){
    in_iquery = .iquerySubsetRedimensionInput(array.in, attributes, min=b$min, max=b$max,  
                                             index.start = index.start, index.end = index.end,
                                             chunk = chunk, overlap)
    proc_iquery = paste("r_exec(",in_iquery,",'output_attrs=",8,"',","'expr=",expr,"')",sep="")
    proc_iquery = .iqueryRedimensionOutput(array.in = proc_iquery, schema = outputschema)
    paste0("insert(",proc_iquery,",",array.out,")")
  })
  res
}

#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an afl query to redimension 
#' a subset of the input array for the TWDTW analysis. 
#' 
#' @inheritParams iqueryCreateRankArray
#' 
#' @seealso \link[dtwSat]{twdtwApplyiquery} 
#' 
#' @examples 
#' ###
#' 
#' @export
iquerySellectMatches = function(array.in, min, max, k=1,
                                chunk, overlap=c(0,0,0,0)){
  
  scidbarray.in = scidb(array.in)
  index.start = as.numeric(scidb_coordinate_start(scidbarray.in))
  index.end = as.numeric(scidb_coordinate_end(scidbarray.in))
  
  time = time[ index.start[3]<=time & time<=index.end[3] ]
  if(length(time)<1)
    stop("time indices out of bounds")
  if(k < index.start[5] & index.end[5] < k )
    stop("Index k out of bounds")
  
  if(missing(min)) min = index.start
  if(missing(max)) max = index.end
  
  col_ids = paste(index.start[1], index.end[1],sep=":")
  row_ids = paste(index.start[2], index.end[2],sep=":")
  time_ids = paste(index.start[3], index.end[3],sep=":")
  class_ids = paste(index.start[4], index.end[4],sep=":")
  
  outputschema = paste("<distance:double>",
                       " [col_id=",col_ids,",",chunk[1],",",overlap[1],
                       ",row_id=",row_ids,",",chunk[2],",",overlap[2],
                       ",time_id=",time_ids,",",index.end[3],",",overlap[3],
                       ",class_id=",time_ids,",",index.end[4],",",overlap[4],"]", sep="")
  
  i_min = c(min[1:4], k)
  i_max = c(max[1:4], k)
  res = paste0("between(",array.in,",",paste(i_min,collapse = ","),",",paste(i_max,collapse = ","),")")
  res = paste0("redimension(",res,",",outputschema,")")
  res
}

#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an afl query to redimension 
#' a subset of the input array for the TWDTW analysis. 
#' 
#' @inheritParams iqueryCreateRankArray
#' 
#' @examples 
#' ###
#' 
#' @export
iqueryBuild2DArray = function(array.in, min, max, t=1,
                                chunk, overlap=c(0,0,0,0)){
  
  scidbarray.in = scidb(array.in)
  index.start = as.numeric(scidb_coordinate_start(scidbarray.in))
  index.end = as.numeric(scidb_coordinate_end(scidbarray.in))
  
  if(t < index.start[3] & index.end[3] < t )
    stop("Time index out of bounds")
  
  if(missing(min)) min = index.start
  if(missing(max)) max = index.end
  
  col_ids = paste(index.start[1], index.end[1],sep=":")
  row_ids = paste(index.start[2], index.end[2],sep=":")
  
  outputschema = paste("<class:uint16,distance:double>",
                       " [col_id=",col_ids,",",chunk[1],",",overlap[1],
                       ",row_id=",row_ids,",",chunk[2],",",overlap[2],"]", sep="")
  
  i_min = c(min[1:2], t, min[4])
  i_max = c(max[1:2], t, max[4])
  res = paste0("between(",array.in,",",paste(i_min,collapse = ","),",",paste(i_max,collapse = ","),")")
  res = paste0("apply(",res,", class, uint16(class_id))") 
  res = paste0("redimension(",res,",",outputschema,")")
  res
}

#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function builds an afl query to redimension 
#' a subset of the input array for the TWDTW analysis. 
#' 
#' @inheritParams iqueryTWDTWBlock
#' 
#' @seealso \link[dtwSat]{twdtwApplyiquery} 
#' 
#' @examples 
#' 
#' 
#' @export
iqueryRankBlock = function(array.in, array.out, grid = NULL, pixelsize = NULL, 
                           blocksize, min, max, chunk, overlap=c(0,0,0,0)){
  
  scidbarray.in = scidb(array.in)
  scidbarray.out = scidb(array.out)
  index.start = as.numeric(scidb_coordinate_start(scidbarray.in))
  index.end = as.numeric(scidb_coordinate_end(scidbarray.in))
  outputschema = schema(scidbarray.out)
  
  if(missing(min)) min = index.start
  if(missing(max)) max = index.end
  if(length(min)<3) min[3] = index.start[3]
  if(length(max)<3) max[3] = index.end[3]
  if(length(min)<4) min[4] = index.start[4]
  if(length(max)<4) max[4] = index.end[4]
  
  col_ids = paste(index.start[1], index.end[1],sep=":")
  row_ids = paste(index.start[2], index.end[2],sep=":")
  time_ids = paste(index.start[3], index.end[3],sep=":")
  class_ids = paste(index.start[4], index.end[4],sep=":")
  
  if(is.null(grid)){
    blocks = .blockSizeFromColRow(min, max, blocksize)
  } else {
    blocks = .blockSizeFromPoly(min, max, grid, pixelsize)
  }
  res = lapply(blocks, function(b){
    in_iquery = .iqueryRankTWDTWDistance(array.in, b$min, b$max, index.start, index.end, chunk, overlap)
    str_iquery = paste0("redimension(",in_iquery,",",outputschema,")")
    paste0("insert(",str_iquery,",",array.out,")")
  })
  res
}


#' @title iquery builder    
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function runs the iquery string on SciDB. 
#' 
#' @param str_iquery A list of iquery strings to be executed in SciBD. 
#' @param logfile The path for the log file. 
#' 
#' @seealso \link[dtwSat]{iqueryBlock} 
#' 
#' @examples 
#' ###
#' 
#' @export
iqueryApply = function(str_iquery, logfile=NULL){
  
  if(is.null(logfile)) logfile = paste0("./log_",format(Sys.time(), "%Y%b%d_%H-%M-%S"))
  log_con = file(logfile, open="a")
  
  query_time = system.time(
    lapply(seq_along(str_iquery), function(i){
      capture.output(Sys.time(), file = log_con)
      cat("Processing block: ", i,"/",length(str_iquery),"\n", file = log_con)
      cat(str_iquery[[i]], "\n\n", file = log_con)
      iquery(str_iquery[[i]], return = FALSE)
    })
  )
  cat("Processing finished at: ", file = log_con)
  capture.output(Sys.time(), file = log_con)
  cat("Processing time in (h): \n", file = log_con)
  capture.output(query_time/60/60, file = log_con)
  close(log_con)
  query_time
}


.blockSizeFromPoly = function(min, max, grid, pixelsize){
  prj_str = projection(grid)
  res = lapply(grid@polygons, function(pol){
    coordinates = data.frame(t(bbox(pol)))
    modis_grid = modisColRowFromLongLat(coordinates, pixelsize = pixelsize, proj4string = prj_str)
    data.frame(min = c(modis_grid$col[1],modis_grid$row[2],min[-c(1,2)]), max=c(modis_grid$col[2],modis_grid$row[1], max[-c(1,2)]))
  })
  res
}

.blockSizeFromColRow = function(min, max, blocksize){
  cols = seq(min[1], max[1], blocksize[1])
  rows = seq(min[2], max[2], blocksize[2])
  k = 1
  res = list()
  for(i in rows)
    for(j in cols){
      res[[k]] = data.frame(min = c(j,i,min[-c(1,2)]), max=c(j+blocksize[1]-1, i+blocksize[1]-1, max[-c(1,2)]))
      if(res[[k]]$max[2] >= max[2])
        res[[k]]$max[2] = max[2]
      if(res[[k]]$max[1] >= max[1])
        res[[k]]$max[1] = max[1]
      k = k + 1
    }
  res
}

.iquerySubsetArray = function(array, attributes, min, max){
  res = paste("between(",array,",",paste(min,collapse = ","),",",paste(max,collapse = ","),")", sep="")
  attributes = c("col_id", "row_id", "time_id", attributes)
  datasets2double = paste(paste0("d", attributes), ", double(",attributes,")", sep="", collapse = ",")
  res = paste("apply(",res,",",datasets2double,")", sep="")
  res = paste("project(",res,",",paste0("d", c(attributes), collapse = ","),")", sep="")
  res
}

.iqueryRedimensionOutput = function(array.in, schema){
  r_exec_output = c("expr_value_0", "expr_value_1", "expr_value_2", "expr_value_3", "expr_value_4", "expr_value_5", "expr_value_6", "expr_value_7")
  names(r_exec_output) = c("col_id", "row_id", "time_id", "class_id", "ts_id", "year", "cdoy", "distance")
  dimensions = c("col_id", "row_id", "time_id", "class_id")
  attributes = c("ts_id", "year", "cdoy", "distance")
  dim2int64 = paste(dimensions, ", int64(",r_exec_output[dimensions],")", sep="", collapse = ",")
  attr2double = paste(attributes, ", double(",r_exec_output[attributes],")", sep="", collapse = ",")
  res = paste("apply(",array.in,",",paste(dim2int64,attr2double,sep=","),")", sep="")
  proj_attr = paste(c(dimensions, attributes), sep="", collapse = ",")
  res = paste("project(",res,",",proj_attr,")", sep="")
  res = paste("redimension(",res,",",schema,")", sep="")
  res 
}


.iquerySubsetRedimensionInput = function(array, attributes, min, max, index.start, 
                                         index.end, chunk, overlap){
  
  col_ids = paste(index.start[1], index.end[1],sep=":")
  row_ids = paste(index.start[2], index.end[2],sep=":")
  time_ids = paste(index.start[3], index.end[3],sep=":")
  if(length(chunk)<3) chunk[3] = index.end[3]
  
  res = .iquerySubsetArray(array, attributes, min, max)
  
  attributes = c("col_id", "row_id", "time_id", attributes)
  datasets2double = paste(paste0("d", attributes), ":double", sep="", collapse = ",")
  
  res = paste("redimension(",res,", <",datasets2double,"> [col_id=",col_ids,",",chunk[1],",",
              overlap[1],",row_id=",row_ids,",",chunk[2],",",overlap[2],",time_id=",time_ids,",",
              chunk[3],",",overlap[3],"])",sep="")
  res
}

.iqueryRankTWDTWDistance = function(array, min, max, index.start, index.end, chunk, overlap){
  
  col_ids = paste(index.start[1], index.end[1],sep=":")
  row_ids = paste(index.start[2], index.end[2],sep=":")
  time_ids = paste(index.start[3], index.end[3],sep=":")
  class_ids = paste(index.start[4], index.end[4],sep=":")
  
  schema_in =  paste("<distance:double>",
                     " [col_id=",col_ids,",",chunk[1],",",0,
                     ",row_id=",row_ids,",",chunk[2],",",0,
                     ",time_id=",time_ids,",",chunk[3],",",overlap[3],
                     ",class_id=",class_ids,",",chunk[4],",",overlap[4],"]", sep="")
  
  res = paste0("redimension(",array,",",schema_in,")")
  res = paste0("between(",array,",",paste(min,collapse = ","),",",paste(max,collapse = ","),")")
  res = paste0("rank(",res, ",distance,time_id,col_id,row_id)")
  res = paste0("apply(",res,",rank_id,int64(distance_rank))")
  res
}

