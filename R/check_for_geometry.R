check_for_geometry <- function(in_sf){
  
  if(!"geometry" %in% names(in_sf)){
    stop("Must contain column: geometry")
  }
}
