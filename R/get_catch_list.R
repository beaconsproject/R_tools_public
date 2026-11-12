### get_catch_list ###
#
#' Get unique sets of catchment values from a table.
#' 
#' Designed for use with wide format tables such as those output by \code{get_upstream_catchments} where each 
#' column name is a unique protected area id, and each row is a catchment id. Catchment values from all columns 
#' included in conservation_areas_ids are combined and returned after removing duplicates.
#' 
#' @param conservation_areas_ids A single column name in input_table, or a vector of multiple column names.
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
get_catch_list <- function(conservation_areas_ids, input_table){
  
  # check nets are in input_table
  if(!all(conservation_areas_ids %in% colnames(input_table))){
    stop(paste0("names do not appear as colnames in the input table: ", conservation_areas_ids[!conservation_areas_ids %in% colnames(input_table)]))
  }
  
  # get unique catchments vector from input_table using vector of colnames
  net_catchments <- unique(unlist(lapply(conservation_areas_ids, function(x){
    input_table[[x]]
  })))
  
  net_catchments[!is.na(net_catchments)]
}
