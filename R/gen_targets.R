### gen_targets ###
#
#' Generate target tables for representation analysis.
#'
#'Takes a reference area polygon and a classified raster for which representation will be evaluated. Calculates the proportion of each class in the reference area 
#'and multiples by reserve size to get target in km2.
#'
#' @param reference_sf sf object of the reference area we are aiming to represent
#' @param representation_raster Raster object of the representation layer classified into categorical classes
#' @param class_values A vector of classes in representation_raster to generate targets for. Defaults to all classes in the representation_raster.
#' @param reserve_size The area in km2 that targets will sum to. Generally the approximate reserve size of reserves being evaluated.
#'
#' @return A tibble with columns: 
#'\itemize{
#'  \item{class_value: the list of class_values}
#'  \item{area_km2: the area of each class_value in the reference_sf}
#'  \item{reserve_size: the provided reserve_size}
#'  \item{class_proportion: area_km2/sum(area_km2)}
#'  \item{target_km2: class_proportion * reserve_area}
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' gen_targets(ref_poly, led_sample, 5000)
gen_targets <- function(reference_sf, representation_raster, reserve_size, class_values = c()){
  
  stopifnot(sf::st_crs(reference_sf) == sf::st_crs(representation_raster))
  
  cellsize <- raster::res(representation_raster)[1] / 1000
  cellarea <- cellsize * cellsize
  
  # if more than one feature in the reference object, dissolve into single geometry
  if(nrow(reference_sf) > 1){
    reference_sf <- reference_sf %>%
      dplyr::summarise(geometry = sf::st_union(.data$geometry))
  }
  
  x <- exactextractr::exact_extract(representation_raster, reference_sf, progress = FALSE)[[1]]
  
  # if no class_values provided, use all values in extracted raster
  if(length(class_values) == 0){
    class_values <- unique(x$value)
    class_values <- class_values[!is.na(class_values)]
  }
  
  # calc sum of each class_value
  df <- x %>%
    dplyr::filter(.data$value %in% class_values) %>%
    dplyr::mutate(area = .data$coverage_fraction * cellarea) %>%
    dplyr::group_by(.data$value) %>%
    dplyr::summarize(area_km2 = sum(.data$area))
  
  names(df) <- c("class_value", "area_km2")
  
  # if any class_values not in df, add them with area of zero
  for(i in class_values){
    if(!i %in% df$class_value){
      df <- rbind(df, data.frame(class_value = i, area_km2 = 0))
    }
  }
  
  # order by class_value
  df <- df[order(df$class_value),]
  
  # add reserve size and calc targets
  df$reserve_size <- reserve_size
  df$class_proportion <- df$area_km2 / sum(df$area_km2)
  df$target_km2 <- round(df$class_proportion * df$reserve_size, 2)
  
  return(df)
}
