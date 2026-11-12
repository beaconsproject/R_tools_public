# convert list of vectors into df with list element names as colnames and missing values as NAs (i.e. BUILDER style table)
list_to_wide <- function(values_list){
  
  # get out table nrow (i.e. longest list of values)
  tbl_rows <- max(unlist(lapply(values_list, function(x){
    length(x)
  })))
  
  values_list_nas <- lapply(values_list, function(x){
    c(x, rep(NA, tbl_rows - length(x)))
  })
  
  out_tab <- as.data.frame(do.call(cbind, values_list_nas))
  
  return(out_tab)
}
