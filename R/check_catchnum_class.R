#check_catchnum <- function(catchments_sf){
  
  # check CATCHNUM exists
#  if(!"CATCHNUM" %in% names(catchments_sf)){
#    stop("Catchments must contain column 'CATCHNUM'")
#  }
#}

check_catchnum_class <- function(catchments_sf, builder_table){
  col_classes <- sapply(colnames(builder_table), function(x) class(builder_table[[x]]))
  if(!all(col_classes == class(catchments_sf$CATCHNUM))){
    warning(paste0("Table column classes do not match class(catchments_sf$CATCHNUM) for columns: ", paste0(colnames(builder_table)[col_classes != class(catchments_sf$CATCHNUM)], collapse=", ")))
  }
}
