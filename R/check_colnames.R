check_colnames <- function(x, x_name, cols){
  for(col in cols){
    if(!col %in% colnames(x)){
      stop(paste0("Column '", col, "' not in table '", x_name, "'"))
    }
  }
}
