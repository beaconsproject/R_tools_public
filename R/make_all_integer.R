make_all_integer <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.integer(x[[col]])){
        warning(paste0("is.integer(", col, ") == FALSE; converting to integer"))
        x[[col]] <- as.integer(x[[col]])
      }
    }
  }
  return(x)
}
