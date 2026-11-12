### raster_mean ###
#
#' Calculate the mean of a raster layer inside a set of input polygons.
#'
#' Clips the raster layer to each conservation area polygon and calculates the mean of the raster values.
#'
#' @param conservation_areas_sf sf object of conservation areas in which to sum raster values.
#' @param raster_layer Raster object that will be clipped to the conservation areas, with crs matching conservation_areas_sf.
#'
#' @return A vector of mean values matching the order of the input \code{conservation_areas_sf}.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network")
#' 
#' raster_mean(conservation_areas, led_sample)

raster_mean <- function(conservation_areas_sf, raster_layer){
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(raster_layer))
  
  # check geometry column is present in sf objects
  check_for_geometry(conservation_areas_sf)
  
  # set up output vector
  result_vector <- c()
  
  # make iterating column for subsetting
  conservation_areas_sf$id <- as.character(1:nrow(conservation_areas_sf))
  
  # split networks into blocks of 10 for processing
  id_list <- conservation_areas_sf$id
  id_list_grouped <- split(id_list, ceiling(seq_along(id_list)/10))
  
  for(id_list_g in id_list_grouped){
    
    conservation_areas_sf_id <- conservation_areas_sf[conservation_areas_sf$id %in% id_list_g,] # subset dissolved networks by the block of networks
    x <- exactextractr::exact_extract(raster_layer, conservation_areas_sf_id, progress = FALSE) # extract
    
    names(x) <- id_list_g # name the list elements by their associated netname
    
    for(id_i in id_list_g){
      
      # for each reserve in the block, extract values and calculate mean
      vals <- x[[id_i]] %>%
        dplyr::filter(!is.na(.data$value)) %>%
        dplyr::filter(.data$coverage_fraction > 0.5) %>% # only keep values from cells with at least half their area in the polygon
        dplyr::pull(.data$value)
      
      result <- mean(vals, na.rm = TRUE)
      
      # add result to return vector
      result_vector <- c(result_vector, round(result, 4))
    }
  }
  return(result_vector)
}
