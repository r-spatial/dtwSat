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
                SM = matrix(as.integer(step.matrix), nrow(step.matrix), ncol(step.matrix)),
                N  = as.integer(n),
                M  = as.integer(m),
                NS = as.integer(nrow(step.matrix)),
                PACKAGE="dtwSat")
  }else{#Take a coffee
        lm = matrix(NA, nrow=n, ncol=m)
        lm[1,] = 0
        
        nsteps = dim(step.matrix)[1]
        lm[1,1] = cm[1,1]
        sm = matrix(NA, nrow=n, ncol=m)

        dir = step.matrix
        ns = attr(dir,"npat")
        trash = lapply(1:m, function(j){
          trash = lapply(1:n, function(i){
            if(!is.na(lm[i,j]))
              return(NULL)
            
            clist = numeric(ns)+NA
            for(k in 1:nsteps){
              p = dir[k,1]
              I = i-dir[k,2]
              J = j-dir[k,3]
              if(I>=1 && J>=1) {
                w = dir[k,4]
                if(w == -1) {
                  clist[p] = lm[I,J]
                }else{
                  clist[p] = clist[p]+w*cm[I,J]
                }
              }
            }
            minc = which.min(clist)
            if(length(minc) > 0){
              lm[i,j] <<- clist[minc]
              sm[i,j] <<- minc
            }
          })
        })
        
        out = list(CM = lm, DM = sm)
   }
   res = list()
   res$costMatrix = out$CM[-1,]
   res$directionMatrix = out$DM[-1,]
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
    }else{# Take a coffee
        res = lapply(jmin, function(j){
              i = n
              dir = step.matrix
              ns = attr(dir,"npat")
              
              nullrows = dir[,2]==0 & dir[,3]==0
              tmp = dir[!nullrows,,drop=FALSE]
              
              steps.list = lapply(1:ns, function(k){
                sbs = tmp[,1]==k  
                spl = tmp[sbs,-1,drop=FALSE]
                nr = nrow(spl)
                spl[nr:1,,drop=FALSE]
              })
              
              I = c(i)
              J = c(j)
              
              repeat{
                    if(i==1)
                      break     
                    s = dm[i,j]
                    if(is.na(s))
                      break
                    
                    steps = steps.list[[s]]
                    ns = nrow(steps)
                    
                    trash = lapply(1:ns, function(k){
                      if(i-steps[k,1] > 0){
                        I <<- c(i-steps[k,1],I)
                        J <<- c(j-steps[k,2],J)
                      }   
                      NULL
                    })
                    
                    i = I[1]
                    j = J[1]
              }
              out = list(index1 = I, index2 = J)
              out
          })
    }
    res
}





