make_all_numeric <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.numeric(x[[col]])){
        warning(paste0("is.numeric(", col, ") == FALSE; converting to numeric"))
        x[[col]] <- as.numeric(x[[col]])
      }
    }
  }
  return(x)
}
