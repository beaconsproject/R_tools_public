### group_conservation_areas ###
#
#' Group conservation areas based on overlap using a grid.
#' 
#' Groups are assigned based on polygon centroid membership within 
#' an intersecting grid. Centroids in the same grid cell are assigned the same group id.
#'
#' @param conservation_areas_sf sf object of conservation areas.
#' @param grid_size Numeric size of grid cells to create in units matching the crs of \code{conservation_areas_sf}.
#'
#' @return Vector of group ids matching the input polygons.
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
#' group_conservation_areas(reserves, 10000)
group_conservation_areas <- function(conservation_areas_sf, grid_size){
  
  sf::st_agr(conservation_areas_sf) = "constant"
  
  # Make centroids - or pointOnPoly
  points_sf <- conservation_areas_sf %>%
    sf::st_centroid()
  
  # Make grid covering points
  grid_sf <- sf::st_make_grid(points_sf, cellsize = c(grid_size, grid_size), what = 'polygons')
  grid_sf <- sf::st_sf(geometry = grid_sf, data.frame('grid_id' = 1:length(grid_sf)))
  
  # remove grid_id if it already exists
  if("grid_id" %in% colnames(points_sf)){
    points_sf <- points_sf %>%
      dplyr::select(-all_of("grid_id"))
  }
  
  # spatial join fishnet grid id to conservation areas
  points_sf <- points_sf %>%
    sf::st_join(grid_sf)
  
  # convert grid ids to ordered vector starting at 1
  mappings <- data.frame(grid_id = unique(points_sf$grid_id), new_val = 1:length(unique(points_sf$grid_id)))
  
  points_sf <- points_sf %>%
    dplyr::left_join(mappings, by = "grid_id")
  
  return(points_sf$new_val)
}