
#' @title Splits the total number of processes
#' 
#' @description A recursive function to splits the total 
#' number of processes according to the threads size.
#' 
#' @param n An integer for the total number of processes. 
#' @param thread.size An integer for the thread size.
#' @docType methods
#' @export
computeThreadsSize = function(n, thread.size){
  if( n < 1.5*thread.size || thread.size < 1)
    return(n)
  return( c(thread.size, computeThreadsSize(n-thread.size, thread.size)) )
}


#' @title Computes a list of GridTopology objects
#' 
#' @description This function computes a list of GridTopology 
#' objects with the cells to be processed in each thread.  
#' 
#' @param cellcentre.offset see ?GridTopology.
#' @param cellsize see ?GridTopology.
#' @param cells.dim see ?GridTopology.
#' @param thread.size integer, vector with number of cells in each thread.
#' @param ncores An integer with the number of cores. Default is 1.
#' 
#' @docType methods
#' @export
computeGridTopologyList = function(cellcentre.offset, cellsize,
                                   cells.dim, thread.size, ncores=1)
  {
  if(missing(thread.size))
    thread.size = c(trunc(cells.dim[1]/ncores), trunc(cells.dim[2]/ncores))

  if(length(thread.size)==1)
    thread.size = rep(thread.size, 2)
  
  threadsSizeX = computeThreadsSize(cells.dim[1], thread.size[1])
  threadsSizeY = computeThreadsSize(cells.dim[2], thread.size[2])
  
  res = list()
  k=1
  y = cellcentre.offset[2]
  for(I in threadsSizeY){
    x = cellcentre.offset[1]
    for(J in threadsSizeX){
      res[[k]] = GridTopology(c(x,y), c(cellsize[1],cellsize[2]), c(J, I))
      x = x + J*cellsize[1]
      k = k + 1
    }
    y = y + I*cellsize[2]
  }
  
  res
}


#' @title Recover MODIS dates
#' 
#' @description The function recovers the MODIS image dates. 
#' 
#' @param years An vector of integers with the years.Default 
#' is form 2000 to the system time year format(Sys.time(), "%Y").
#' @param frequency An integer with the frequency in days.
#' @docType methods
#' @export
recoverMODISDates = function(years, frequency=16){
  if(missing(years))
    years = 2000:format(Sys.time(), "%Y")
  
  dates = unlist(lapply(years, function(y){  
    days = seq(from = 0, to = 365, by = frequency)
    as.Date(days, origin = as.Date(paste(y,"-01-01",sep="")))
  }))
  dates = as.Date(dates, origin="1970-01-01")
  return(dates)
}




#' @title Computes DTW distance
#' 
#' @description The function computes the whole possible alignments 
#' for a set of patterns in one time series. 
#' 
#' @param template A zoo object with the template time series. 
#' It must be larger than the query.
#' @param TemporalPatterns.list A list of patterns. Each node of the list 
#' has two columns data.frame with the dates and the value. 
#' @docType methods
#' @export
computeDTWForAllPatterns = function(template, TemporalPatterns.list, ... ){
  res = lapply(seq_along(TemporalPatterns.list), function(j){
    patternName = names(TemporalPatterns.list)[j]
    x = as.numeric(TemporalPatterns.list[[j]][,2])
    tx = as.Date(TemporalPatterns.list[[j]][,1], origin="1970-01-01")
    query = zoo(x, tx)
    out = timeSeriesAnalysis(query, template, theta=THETA, span=SPAN,
                             normalize=NORMALIZE, satStat=SATSTAT)
    return(out)
  })
  names(res) = names(TemporalPatterns.list)
  return(res)
}


#' @title Build SciDB query
#' 
#' @description The function builds the SciDB query for time 
#' series analysis. 
#' @param INPUTARRAY
#' @param OUTPUTSCHEMA
#' @param PATTERNNAMES
#' @param XMIN
#' @param XMAX
#' @param YMIN
#' @param YMAX
#' @param COLIDS
#' @param ROWIDS
#' @param FROM
#' @param TO
#' @param THETA = 1.0
#' @param SPAN = 0.3
#' @param OVERLAPPING = 0.15
#' @param THRESHOLD = 1.0
#' @param NORMALIZE = TRUE
#' @param SATSTAT = FALSE 
#' @param NCORES = 1
#' @param LIBRARYPATH
#' @param EXPPRPATH
#' @docType methods
#' @export
buildSciDBDTWQuery = function(INPUTARRAY, OUTPUTSCHEMA, PATTERNNAMES, 
                              XMIN, XMAX,
                              YMIN, YMAX,
                              COLIDS, ROWIDS,
                              FROM, TO,
                              THETA = 1.0, SPAN = 0.3,
                              OVERLAPPING = 0.15, THRESHOLD = 1.0,
                              NORMALIZE = TRUE, SATSTAT = FALSE, 
                              NCORES = 1,
                              LIBRARYPATH, EXPPRPATH){

rOut = 0:(length(PATTERNNAMES)+2)
attrRename = gsub(" ", ",",paste(paste("expr_value", rOut, sep="_"), c("col", "row", "year", PATTERNNAMES), collapse = ","))

attrRedimension = paste(c("col_id", "row_id", "year_id", PATTERNNAMES), collapse = ",")

tmin = as.integer(as.Date(FROM, origin="1970-01-01")) - as.integer(as.Date("2000-02-18", origin="1970-01-01") )
tmax = trunc( (as.integer(as.Date(TO, origin="1970-01-01")) - as.integer(as.Date("2000-02-18", origin="1970-01-01")) ) / 23 ) + 1

if(tmin < 0)
  tmin = 0

res = paste(" redimension(
                        project(
                              apply(
                                    attribute_rename(
                                              r_exec(
                                                    project(
                                                              apply(
                                                                    redimension(
                                                                                between(",INPUTARRAY,",",XMIN,",",YMIN,",",tmin,",",XMAX,",",YMAX,",",tmax,"),
                                                                                <evi:int16>[col_id=",COLIDS,",32,0,row_id=",ROWIDS,",32,0,time_id=0:9200,",tmax+1,",0]
                                                                    ),
                                                                    devi, double(evi), dcol, double(col_id), drow, double(row_id), dtime, double(time_id)
                                                            ), 
                                                              devi, dcol, drow, dtime
                                                    ),
                                                    'output_attrs=",length(rOut),"',
                                                    'expr=
                                                          LIBRARYPATH=\"",LIBRARYPATH,"\"
                                                          FROM=\"",FROM,"\"
                                                          TO=\"",TO,"\"
                                                          THETA=",THETA,"
                                                          SPAN=",SPAN,"
                                                          OVERLAPPING=",OVERLAPPING,"
                                                          THRESHOLD=",THRESHOLD,"
                                                          NORMALIZE=",NORMALIZE,"
                                                          SATSTAT=",SATSTAT,"
                                                          NCORES=",NCORES,"
                                                          source(\"",EXPPRPATH,"\")
                                                          res=list(59787,49064,2002,4,5,6,7,8,9,10,11,12,13,14)
                                                          res'
                                             ),
                                             ",attrRename,"
                                     ),
                                     col_id,  int64(col),
                                     row_id,  int64(row),
                                     year_id, int64(year)
                              ),",
                              attrRedimension,"
                        ),",
                        OUTPUTSCHEMA,"
              )", sep="")

return(res)

}

