### group_conservation_area_networks ###
#
#' Group networks based on overlap of their component conservation areas. 
#' 
#' Individual conservation areas are first grouped using [group_conservation_areas()]. Networks are then grouped 
#' based on the combination of their component conservation areas. Two networks with component conservation areas in matching groups will be grouped together.
#' 
#' e.g. for two networks PA1__PA2 and PB3__PB4, if PB1 and PB3 fall in the same grid cell, and PB2 and PB4 fall
#' in the same grid cell, networks will receive the same unique group id.
#'
#' @param network_list Vector of network names to group. Must be in the format x1__x2__x3 where x1, x2 and x3 are 
#' values in the \code{conservation_areas_sf} \code{conservation_areas_id} column. Any number of x can be used (e.g. x1__x2, 
#' x1__x2__x3__x4 etc.).
#' @param conservation_areas_sf sf object of conservation areas.
#' @param conservation_areas_id unique identifier column in \code{conservation_areas_sf}.
#' @param grid_size Numeric size of grid cells to create in units matching the crs of \code{conservation_areas_sf}.
#'
#' @return Vector of group ids matching \code{network_list}.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = colnames(builder_table_sample))
#' network_names <- gen_network_names(conservation_areas$network, 2)
#' networks_sf <- build_network_polygons(conservation_areas, network_names)
#' 
#' group_conservation_area_networks(networks_sf$network, conservation_areas, "network", 10000)

group_conservation_area_networks <- function(network_list, conservation_areas_sf, conservation_areas_id, grid_size){
  
  # get all conservation areas names, check they are all in conservation_areas_sf, subset conservation_areas_sf
  ca_names <- unique(unlist(sep_network_names(network_list)))
  
  if(all(ca_names %in% conservation_areas_sf[[conservation_areas_id]])){
    
    conservation_areas_sf <- conservation_areas_sf %>%
      dplyr::filter(.data[[conservation_areas_id]] %in% ca_names)
  } else{
    stop(paste0("Conservation area missing from conservation_areas_sf: ", ca_names[!ca_names %in% conservation_areas_sf[[conservation_areas_id]]]))
  }
  
  # group
  conservation_areas_sf$grid_id <- group_conservation_areas(conservation_areas_sf, grid_size)
  
  # make tibble of network names
  tb <- dplyr::tibble(network_list = network_list)
  
  # calculate group_id combinations and add to tibble
  tb$grid_id <- sapply(network_list, function(x){
    nets <- sep_network_names(x)
    paste0(sort(sapply(nets, function(y){conservation_areas_sf$grid_id[conservation_areas_sf[[conservation_areas_id]] == y]})), collapse = "_")
  })
  
  # make table of mappings for unique grid code combinations
  mappings <- data.frame(grid_id = unique(tb$grid_id), new_val = 1:length(unique(tb$grid_id)))
  
  # join output grid ids to tb
  tb <- tb %>%
    dplyr::left_join(mappings, by = "grid_id")
  
  return(tb$new_val)
}