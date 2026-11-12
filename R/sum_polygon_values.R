### sum_polygon_values ###
#
#' Sum area of a polygon layer inside a set of input polygons.
#'
#' Sums the area of a polygon dataset (or of each group within the polygon dataset) within each conservation area polygon.
#' 
#' Sums are returned in a tibble which can be joined to \code{conservation_areas_sf} using a \code{dplyr::left_join()}
#'  and \code{dplyr::pivot_wider()}.
#'  
#' Projections of \code{conservation_areas_sf} and \code{sum_polygon} must match and are assumed to have units in meters.
#'
#' @param conservation_areas_sf sf object of conservation areas in which to sum raster values.
#' @param conservation_areas_id String representing unique id column in \code{conservation_areas_sf}.
#' @param sum_polygon sf object of one or more polygons to sum.
#' @param sum_polygon_groups Column name in \code{sum_polygon} holding grouping variables. If provided, area of each group is summed, otherwise
#' full area of \code{sum_polygon} is summed. 
#' @param class_vals Vector of class values in \code{sum_polygon_groups} to include in the output. Defaults to all non-NA values in \code{sum_polygon_groups}.
#' Ignored if \code{sum_polygon_groups} is not provided.
#' @param fill_zeros Should missing values be included with an area of zero? Default is TRUE.
#'
#' @return A tibble listing area in km2 for \code{sum_polygon} or \code{sum_polygon_groups} in each conservation area.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library(dplyr)
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network")
#' 
#' # Make a habitat polygon with 2 groups using catchments  
#' habitat <- catchments_sample %>%
#'              filter(CATCHNUM %in% c(191396, 191649)) %>%
#'              mutate(group = c(1,2)) %>%
#'              select(group)
#'              
#' sum_polygon_values(conservation_areas, "network", habitat)
#' sum_polygon_values(conservation_areas, "network", habitat, "group")
#' sum_polygon_values(conservation_areas, "network", habitat, "group", c(1))
#' sum_polygon_values(conservation_areas, "network", habitat, "group", c(2))
#' sum_polygon_values(conservation_areas, "network", habitat, "group", c(2), fill_zeros = FALSE)
#' sum_polygon_values(conservation_areas, "network", habitat, "group", fill_zeros = FALSE)
#' 
#' # To join results back to conservation_areas
#' library(tidyr)
#' library(dplyr)
#' conservation_areas %>%
#'   left_join(
#'     pivot_wider(
#'       sum_polygon_values(conservation_areas, "network", habitat, "group"), 
#'       names_from = group, values_from = area_km2)
#'    )
sum_polygon_values <- function(conservation_areas_sf, conservation_areas_id, sum_polygon, sum_polygon_groups = NULL, class_vals = c(), fill_zeros = TRUE){
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(sum_polygon))
  
  # check geometry column is present in sf objects
  check_for_geometry(conservation_areas_sf)
  check_for_geometry(sum_polygon)
  
  sf::st_agr(conservation_areas_sf) = "constant"
  sf::st_agr(sum_polygon) = "constant"
  
  if(!is.null(sum_polygon_groups)){
    
    # check sum_polygon_groups is in sum_polygon
    if(!sum_polygon_groups %in% names(sum_polygon)){
      stop(paste0(sum_polygon_groups, " is not a column in sum_polygon"))
    }
    
    # Grouped code
    out_df <- dplyr::tibble(a = as.character(), b = as.character(), area_km2 = as.numeric()) # make table template
    names(out_df)[names(out_df) %in% c("a","b")] <- c(conservation_areas_id, sum_polygon_groups)
    
    sum_polygon_dslv <- sum_polygon %>%
      dplyr::group_by(.data[[sum_polygon_groups]]) %>%
      sf::st_union() # make sure sum_polygon is a single feature per group
    
    all_groups <- unique(sum_polygon[[sum_polygon_groups]]) # get list of all groups
    
    for(reserve in conservation_areas_sf[[conservation_areas_id]]){ # loop conservation_areas and intersect grouped sum_polygon
      
      # need to split this operation into two so the grouped sum_polygons can be made "constant". Otherwise we get the warning msg.
      df <- sum_polygon %>%
        dplyr::group_by(.data[[sum_polygon_groups]])
      sf::st_agr(df) = "constant"
      
      df <- df %>%
        sf::st_intersection(conservation_areas_sf[conservation_areas_sf[[conservation_areas_id]] == reserve,]) %>%
        dplyr::mutate(area_km2 = round(as.numeric(sf::st_area(.data$geometry)) / 1000000, 3)) %>% # assumes units of m
        sf::st_drop_geometry() %>%
        dplyr::select(.data[[conservation_areas_id]], .data[[sum_polygon_groups]], .data$area_km2) %>%
        dplyr::filter(.data$area_km2 > 0) %>% # intersection includes touching polygons that don't actually overlap. Don't want to return these
        dplyr::as_tibble()
      
      if(fill_zeros){
        # fill in missing values with 0
        for(group_id in all_groups[!all_groups %in% df[[sum_polygon_groups]]]){
          add_row <- dplyr::tibble(a = reserve, b = group_id, area_km2 = 0)
          names(add_row)[names(add_row) %in% c("a","b")] <- c(conservation_areas_id, sum_polygon_groups)
          df <- rbind(df, add_row)
        }
      }
      out_df <- rbind(out_df, df)
      
      # filter by class_vals if requested
      if(length(class_vals > 0)){
        out_df <- out_df %>%
          dplyr::filter(.data[[sum_polygon_groups]] %in% class_vals)
      }
    }
    
  } else{
    # Non grouped code
    out_df <- dplyr::tibble(a = as.character(), area_km2 = as.numeric())
    names(out_df)[names(out_df) == "a"] <- c(conservation_areas_id)
    
    sum_polygon_dslv <- sum_polygon %>%
      sf::st_union() # make sure sum_polygon is a single feature
    
    for(reserve in conservation_areas_sf[[conservation_areas_id]]){ # loop conservation_areas and intersect sum_polygon
      
      df <- conservation_areas_sf[conservation_areas_sf[[conservation_areas_id]] == reserve,] %>%
        sf::st_intersection(sum_polygon_dslv) %>%
        dplyr::mutate(area_km2 = round(as.numeric(sf::st_area(.data$geometry)) / 1000000, 3)) %>% # assumes units of m
        sf::st_drop_geometry() %>%
        dplyr::select(.data[[conservation_areas_id]], .data$area_km2) %>%
        dplyr::filter(.data$area_km2 > 0) %>% # intersection includes touching polygons that don't actually overlap. Don't want to return these
        dplyr::as_tibble()
      
      # Fill in zeros if missing (i.e. if the polygons don't intersect)
      if(fill_zeros & nrow(df) == 0){
        df <- dplyr::tibble(a = reserve, area_km2 = 0)
        names(df)[names(df) == "a"] <- c(conservation_areas_id)
      }
      out_df <- rbind(out_df, df)
    }
  }
  return(out_df)
}