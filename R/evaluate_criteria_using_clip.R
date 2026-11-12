### evaluate_criteria_using_clip ###
#
#' Evaluate criteria by clipping with a conservation areas or a network.
#'
#' Clips the criteria raster to the conservation area and calculates the proportion of each class within it. If a unique 
#' identifer (conservation_areas_id) is provided, proportions are computed separately for each corresponding polygon. 
#'
#' @param conservation_areas_sf sf object of the reference area we are aiming to represent
#' @param conservation_areas_id Unique id column in conservation_areas_sf containing conservation area names as strings. Can be a unique identifier representing a network. 
#' @param criteria_raster Raster object of the representation layer classified into categorical classes
#' @param class_values A vector of classes in representation_raster to generate targets for. Defaults to all classes in the representation_raster.
#' @param target_size The area in km2 that targets will sum to. Generally the approximate target size of reference area being evaluated. set target Must be set to TRUE

#' @return A tibble with columns: 
#'\itemize{
#'  \item{class_value: the list of class_values}
#'  \item{area_km2: the area of each class_value in the conservation_areas_sf}
#'  \item{class_proportion: area_km2/sum(area_km2)}
#'  
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' evaluate_criteria_using_clip(ref_poly, led_sample)
evaluate_criteria_using_clip <- function(conservation_areas_sf, criteria_raster, conservation_areas_id = NULL, class_values = c(), target_size = NULL){
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(criteria_raster))

  ## Dissolve geometry
  if (!is.null(conservation_areas_id) && conservation_areas_id %in% names(conservation_areas_sf)) {
    conservation_areas_sf <- conservation_areas_sf |>
      dplyr::group_by(.data[[conservation_areas_id]]) |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
  #} else if (nrow(conservation_areas_sf) > 1L) {
  } else {
    conservation_areas_sf <- conservation_areas_sf |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
  }
  
  cellsize <- raster::res(criteria_raster)[1] / 1000
  cellarea <- cellsize * cellsize
  
  x <- exactextractr::exact_extract(criteria_raster, conservation_areas_sf, include_cols = conservation_areas_id, progress = FALSE)
  
  x_df <- purrr::map_dfr(x, ~ .)
  
  # if no class_values provided, use all values in extracted raster
  if(length(class_values) == 0){
    class_values <- unique(x_df$value)
    class_values <- class_values[!is.na(class_values)]
  }
  
  group_vars <- if (!is.null(conservation_areas_id) && conservation_areas_id != "") c(conservation_areas_id, "value") else "value"
  
  # calc sum of each class_value
  df <- x_df %>%
    dplyr::filter(value %in% class_values) %>%
    dplyr::mutate(area = coverage_fraction * cellarea) %>%
    dplyr::group_by_at(group_vars) %>%
    #dplyr::group_by(!!rlang::sym(conservation_areas_id), value) %>%
    dplyr::summarize(area_km2 = sum(area), .groups = "drop")
  
  names(df) <- c(conservation_areas_id, "class_value", "area_km2")
  
  # if any class_values not in df, add them with area of zero
  for(i in class_values){
    if(!i %in% df$class_value){
      df <- rbind(df, data.frame(class_value = i, area_km2 = 0))
    }
  }
  
  # order by class_value and calculate class_proportion
  if (!is.null(conservation_areas_id) && conservation_areas_id != "") {
    df <- df %>%
      dplyr::arrange(!!rlang::sym(conservation_areas_id), class_value) %>%
      dplyr::group_by(!!rlang::sym(conservation_areas_id)) %>%
      dplyr::mutate(class_proportion = area_km2 / sum(area_km2)) %>%
      dplyr::ungroup()
  } else {
    df <- df %>%
      dplyr::arrange(class_value) %>%
      dplyr::mutate(class_proportion = area_km2 / sum(area_km2))
  }
  
  if(!is.null(target_size)){
    df$target_size <- target_size
    df$target_km2 <- round(df$class_proportion * df$target_size, 2)
  }
  
  return(df)
}
