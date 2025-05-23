make_catchnum_integer <- function(catchments_sf){
  
  # check CATCHNUM exists
  if(!"CATCHNUM" %in% names(catchments_sf)){
    stop("Catchments must contain column 'CATCHNUM'")
  }
  
  catchments_sf$CATCHNUM <- as.integer(catchments_sf$CATCHNUM)
  
  return(catchments_sf)
}

check_for_geometry <- function(in_sf){
  
  if(!"geometry" %in% names(in_sf)){
    stop("Must contain column: geometry")
  }
}

logical_to_integer <- function(x){
  
  if(!x %in% c(TRUE, FALSE)){
    stop("input must be TRUE or FALSE")
  }
  
  return(as.integer(x))
}

check_seeds_areatargets <- function(seeds){
  if(!all(seeds$Areatarget > 0)){
    stop("All Areatarget values must be > 0")
  }
}

check_seeds_in_catchments <- function(seeds, catchments_sf) {
  # Make sure CATCHNUM is character before splitting
  seeds_clean <- seeds %>%
    mutate(CATCHNUM = as.character(CATCHNUM)) %>%
    separate_rows(CATCHNUM, sep = ",") %>%
    mutate(CATCHNUM = trimws(CATCHNUM))  # Remove any whitespace
  
  # Convert to same type as in catchments_sf
  catchnum_vals <- catchments_sf$CATCHNUM
  if (is.numeric(catchnum_vals)) {
    seeds_clean$CATCHNUM <- as.numeric(seeds_clean$CATCHNUM)
  }
  
  # Perform check
  if (!all(seeds_clean$CATCHNUM %in% catchnum_vals)) {
    if (any(seeds_clean$CATCHNUM %in% catchnum_vals)) {
      warning("Not all seeds are in catchments_sf")  # some overlap
    } else {
      stop("None of the seeds are in catchments_sf")  # no overlap
    }
  }
}

check_colnames <- function(x, x_name, cols){
  for(col in cols){
    if(!col %in% colnames(x)){
      stop(paste0("Column '", col, "' not in table '", x_name, "'"))
    }
  }
}

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

make_all_character <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.character(x[[col]])){
        warning(paste0("is.character(", col, ") == FALSE; converting to character"))
        x[[col]] <- as.character(x[[col]])
      }
    }
  }
  return(x)
}