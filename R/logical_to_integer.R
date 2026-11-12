logical_to_integer <- function(x){
  
  if(!x %in% c(TRUE, FALSE)){
    stop("input must be TRUE or FALSE")
  }
  
  return(as.integer(x))
}
