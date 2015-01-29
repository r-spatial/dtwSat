
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



