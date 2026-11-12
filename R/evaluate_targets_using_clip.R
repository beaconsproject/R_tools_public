### evaluate_targets_using_clip ###
#
#' Evaluate targets by clipping.
#'
#' Similar to [evaluate_targets_using_catchments()] but instead of pre-calculating the areas of the criteria raster in the catchments,
#' this method clips the criteria raster directly to each provided polygon.
#' Generally used for evaluating conservation areas that did not come from the beaconsbuilder package and do not conform to catchment boundaries.
#' Slower than [evaluate_targets_using_catchments()] for large numbers of polygons.

#' @param conservation_areas_sf sf object containing conservation areas to evaluate
#' @param conservation_areas_id Unique id column in conservation_areas_sf containing conservation area names as strings.
#' @param representation_raster Raster object of the criteria layer that will be evaluated, with crs matching conservation_areas_sf
#' @param target_table Data frame containing columns "class_value", "class_proportion" and "target_km2". i.e. the output from [gen_targets()]. All classes in the target table are 
#'   evaluated. All class_values must match the values in the representation_raster
#' @param conservation_areas_list The conservation_areas_id values to process. Defaults to all conservation_areas_id's in conservation_areas_sf
#'
#' @return A tibble with columns: 
#'\itemize{
#'  \item{\code{class_value}: the list of class_values from the target_table}
#'  \item{\code{area_km2}: area of the class value in the network}
#'  \item{\code{network}: name of the network copied from network_list}
#'  \item{\code{class_proportion}: copied from target_table}
#'  \item{\code{target_km2}: copied from target_table}
#'  \item{\code{prop_target_met}: area_km2/target_km2}
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' target_table <- gen_targets(ref_poly, led_sample, 1600)
#' pas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network",
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' evaluate_targets_using_clip(pas, "network", led_sample, target_table)
evaluate_targets_using_clip <- function(conservation_areas_sf, conservation_areas_id, representation_raster, target_table, conservation_areas_list=c()){
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(representation_raster))
  
  # set conservation_areas_list if not provided
  if(length(conservation_areas_list) == 0){
    conservation_areas_list <- conservation_areas_sf[[conservation_areas_id]]
  }
  
  cell_area <- prod(raster::res(representation_raster)) / 1000000 # convert to area in km2, assumes raster res is in metres
  
  counter <- 1
  for(reserve in conservation_areas_list){
    
    # STEP 1 - sum area of each class in each network
    # subset
    reserve_sf <- conservation_areas_sf[conservation_areas_sf[[conservation_areas_id]] == reserve,]
    
    # clip raster
    x <- exactextractr::exact_extract(representation_raster, reserve_sf, progress = FALSE)[[1]]
    names(x) <- c("class_value", "coverage_fraction")
    rows <- x %>%
      dplyr::filter(.data$class_value %in% target_table$class_value) %>%
      dplyr::mutate(area = .data$coverage_fraction * cell_area) %>%
      dplyr::group_by(.data$class_value) %>%
      dplyr::summarize(area_km2 = sum(.data$area))
    
    missing_class_values <- target_table$class_value[!target_table$class_value %in% rows$class_value]
    
    rows <- rows %>% 
      dplyr::add_row(class_value = missing_class_values, area_km2 = 0) %>%
      dplyr::mutate(network = reserve) %>%
      dplyr::arrange(.data$class_value)
    
    if(counter == 1){
      df <- rows
      counter <- counter + 1
    } else{
      df <- rbind(df, rows)
    }
  }
  
  # STEP 2 - join targets and calculate proportion of target met
  df <- df %>%
    dplyr::left_join(target_table[c("class_value","class_proportion","target_km2")], by = "class_value")
  
  # add proportion of target met
  df$prop_target_met <- round(df$area_km2 / df$target_km2, 2)
  
  return(df)
}
