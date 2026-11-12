### extract_catchments_from_table ###
#
#' Get list of catchments from a table and extract to a new feature.
#' 
#' Takes a list of catchments from a table column using the column name, and extracts the 
#' catchments from the catchments dataset. If multiple column names are provided, the unique combination of 
#' catchments from all columns is returned.
#' 
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM.
#' @param input_table Data frame where column names are polygon names and rows are catchments making up each polygon.
#' e.g. output from [get_upstream_catchments()].
#' @param extract_feature_id String representing the column name to extrat from \code{input_table}. Can also take
#' a vector of names in which case multiple columns will be extracted.
#' @param out_feature_id String representing the output column name holding the polygon unique identifier.
#'   
#' @return A sf object of catchments with unique identifier column \code{out_feature_id}.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' 
#' extract_catchments_from_table(catchments_sample, builder_table_sample, "PB_0001", "network")
#' extract_catchments_from_table(
#'   catchments_sample, builder_table_sample, c("PB_0001", "PB_0002"), "network")
extract_catchments_from_table <- function(catchments_sf, input_table, extract_feature_id, out_feature_id){
  
  check_colnames(catchments_sf, cols = "CATCHNUM") # check for CATCHNUM
  check_catchnum_class(catchments_sf, input_table) # Check catchments match, warning if not
  check_for_geometry(catchments_sf)
  
  # get list of catchments
  catchments_list <- get_catch_list(extract_feature_id, input_table)
    
  # extract
  ex <- catchments_sf %>%
    dplyr::filter(.data$CATCHNUM %in% catchments_list) %>%
    dplyr::mutate(id = paste0(extract_feature_id,collapse="__")) %>%
    dplyr::select(.data$CATCHNUM,.data$id)
    
  # set out name
  names(ex)[names(ex) == "id"] <- out_feature_id
    
  return(ex)
}
