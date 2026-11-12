### get_downstream_catchments ###
#
#' Calculate downstream catchments for a set of input polygons.
#' 
#' Calculates all downstream catchments for each provided protected area polygon and returns as a table, with column names as the unique id 
#' of the protected areas.
#'
#' @param conservation_areas_sf sf object of protected area polygons.
#' @param conservation_areas_id String matching the unique identifier column in \code{conservation_areas_sf}.
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM .
#'
#' @return Tibble where each column name is a unique protected area id, and each row is a catchment making up the 
#' upstream area for that protected area. Blank rows are filled with NA. CATCHNUMs are returned as integers.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' reserves <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
#' get_downstream_catchments(reserves, "network", catchments_sample)
get_downstream_catchments <- function(conservation_areas_sf, conservation_areas_id, catchments_sf){
  
  if(!all(c("ORDER1", "ORDER2", "ORDER3", "BASIN", "CATCHNUM") %in% colnames(catchments_sf))){
    stop("catchments_sf must have attributes: ORDER1, ORDER2, ORDER3, BASIN, CATCHNUM")
  }
  
  # make sure catchnums are integer. Out tables in BUILDER wide format should use integer class for CATCHNUM
  catchments_sf <- make_all_integer(catchments_sf, cols = "CATCHNUM")
  
  # get list of catchnums in each PA
  CAs_catchnums_tab <- catchnums_in_polygon(conservation_areas_sf, conservation_areas_id, catchments_sf)
  
  down_agg_list <- list()
  for(col_id in colnames(CAs_catchnums_tab)){
    
    # get list of catchments
    agg_catchments <- get_catch_list(col_id, CAs_catchnums_tab)
    down_agg <- getAggregationDownstreamCatchments_R(catchments_sf, agg_catchments)
    
    # add to up_agg_list
    down_agg_list[[col_id]] <- down_agg
  }
  
  out_df <- dplyr::as_tibble(list_to_wide(down_agg_list))
  
  return(out_df)
}
