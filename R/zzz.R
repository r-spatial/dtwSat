



.onAttach <- function(lib, pkg)  {
  
  packageStartupMessage(sprintf("Loaded dtwSat v%s. See ?dtwSat for help, citation(\"dtwSat\") for use in publication.\n",
                                utils::packageDescription("dtwSat")$Version ) );
  
  

}

## 
## 
#' @title dtw symmetric step with normalization N
#' 
#' @description A symmetric step with normalization N. 
#' stepPattern object from package dtw, see ?stepPattern.
#' @usage symmetric0 # see ?stepPattern
#' @format symmetric0
#' @export
symmetric0 = dtw:::stepPattern(c(
  1,1,1,-1,
  1,0,0,1,
  2,0,1,-1,
  2,0,0,1,
  3,1,0,-1,
  3,0,0,1
), "N")




