check_catchnum <- function(catchments_sf){
  
  # check CATCHNUM exists
  if(!"CATCHNUM" %in% names(catchments_sf)){
    stop("Catchments must contain column 'CATCHNUM'")
  }
}

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

# remove OID from tables comgin out of BUILDER
remove_oid <- function(in_table){
  if("OID" %in% colnames(in_table)){
    out_table <- in_table %>%
      dplyr::select(-OID)
  } else{
    out_table <- in_table
  }
  return(out_table)
}

# Check for rows
check_for_rows <- function(in_table){
  if(nrow(in_table) == 0){
    stop("input_table has no data")
  }
}

# sep_network_names
#' Split network names into their constituent parts.
#'
#' Takes a vector of network names in the form output by [gen_network_names()] (i.e. \code{"__"} separated),
#' and splits the network names back into their original names.
#' 
#' If a single name is provided, a vector of names is returned. If a vector with length >1 is provided, a list
#' of vectors is returned with each list element named with the network name.
#'
#' @param network_names String of a network name, or a vector of multiple network names. Must include the \code{"__"} separator.
#'
#' @return A vector of names, or a list of vectors.
#' @export
#'
#' @examples
#' sep_network_names("PB_1__PB_2")
#' sep_network_names(c("PB_1__PB_2", "PB_1__PB_3"))
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2)
#' sep_network_names(x)
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 4)
#' sep_network_names(x)

sep_network_names <- function(network_names){
  
  if(length(network_names) == 1){
    out_val <- strsplit(network_names, "__")[[1]]
  }
  
  if(length(network_names) > 1){
    out_val <- lapply(network_names, function(x){
      strsplit(x, "__")[[1]]
    })
    names(out_val) <- network_names
  }
  return(out_val)
}

### get_catch_list ###
#
#' Get unique sets of catchment values from a table.
#' 
#' Designed for use with wide format tables such as those output by \code{get_upstream_catchments} where each 
#' column name is a unique protected area id, and each row is a catchment id. Catchment values from all columns 
#' included in CAs_ids are combined and returned after removing duplicates.
#' 
#' @param pa_ids A single column name in input_table, or a vector of multiple column names.
#' @param input_table Data frame where column names are protected area names and rows are catchments making 
#' up each polygon. e.g. output from [get_upstream_catchments()].
#'   
#' @return A vector of unique catchment values.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' get_catch_list("PB_0001", builder_table_sample)
#' get_catch_list(c("PB_0001", "PB_0003"), builder_table_sample)
get_catch_list <- function(CAs_ids, input_table){
  
  # check nets are in input_table
  if(!all(CAs_ids %in% colnames(input_table))){
    stop(paste0("names do not appear as colnames in the input table: ", CAs_ids[!CAs_ids %in% colnames(input_table)]))
  }
  
  # get unique catchments vector from input_table using vector of colnames
  net_catchments <- unique(unlist(lapply(CAs_ids, function(x){
    input_table[[x]]
  })))
  
  net_catchments[!is.na(net_catchments)]
}