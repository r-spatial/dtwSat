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
#   call Fortran DTW inplementation - 2015-10-27              #
#                                                             #
###############################################################


#' @useDynLib dtwSat computecost
.computecost = function(cm, step.matrix){
  
  cm = rbind(0, cm)
  n = nrow(cm)
  m = ncol(cm)
  
  if(is.loaded("computecost", PACKAGE = "dtwSat", type = "Fortran")){
    out = .Fortran("computecost", 
                   CM = matrix(as.double(cm), n, m),
                   DM = matrix(as.integer(0), n, m),
                   VM = matrix(as.integer(0), n, m),
                   SM = matrix(as.integer(step.matrix), nrow(step.matrix), ncol(step.matrix)),
                   N  = as.integer(n),
                   M  = as.integer(m),
                   NS = as.integer(nrow(step.matrix)),
                   PACKAGE="dtwSat")
  } else {
    stop("Fortran computecost lib is not loaded")
  }
  
  res = list()
  res$costMatrix = out$CM[-1,]
  res$directionMatrix = out$DM[-1,]
  res$startingMatrix = out$VM[-1,]
  res$stepPattern = step.matrix
  res$N = n - 1
  res$M = m
  res
}


#' @useDynLib dtwSat tracepath
.tracepath = function(dm, step.matrix, jmin){

    n = nrow(dm)
    m = ncol(dm)
    if(is.null(jmin))
      jmin = m
    
    if(is.loaded("tracepath", PACKAGE = "dtwSat", type = "Fortran")){
      aloc = length(jmin)*10*n
      paths = .Fortran("tracepath", 
                       DM   = matrix(as.integer(dm), n, m),
                       SM   = matrix(as.integer(step.matrix), nrow(step.matrix), ncol(step.matrix)),
                       JMIN = as.vector(as.integer(jmin)),
                       IND1 = rep(as.integer(0), aloc),
                       IND2 = rep(as.integer(0), aloc),
                       POS  = as.vector(rep(as.integer(0),length(jmin)+1)),
                       N    = as.integer(n),
                       M    = as.integer(m),
                       NS   = as.integer(nrow(step.matrix)),
                       NJ   = as.integer(length(jmin)),
                       AL   = as.integer(aloc),
                       PACKAGE="dtwSat")
      
      res = lapply(seq_along(paths$POS)[-1], function(p){
        I = paths$POS[p]:((paths$POS[p-1])+1)
        list(index1 = paths$IND1[I], index2 = paths$IND2[I])
      })
    }else{
      stop("Fortran tracepath lib is not loaded")
    }
    
    res
}





