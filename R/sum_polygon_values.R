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

  # check required columns
  check_for_geometry(conservation_areas_sf)
  check_for_geometry(sum_polygon)
  missing_cols <- check_colnames(conservation_areas_sf, conservation_areas_id)
  if (length(missing_cols) > 0) stop("Column(s) missing from conservation_areas_sf: ", paste(missing_cols, collapse = ", "))
  if (!is.null(sum_polygon_groups) && length(check_colnames(sum_polygon, sum_polygon_groups)) > 0) {
    stop(sum_polygon_groups, " is not a column in sum_polygon")
  }

  # 1. Merge sum_polygon into one feature per group (or a single feature if no groups)
  if (is.null(sum_polygon_groups)) {
    grp <- ".all"
    sum_polygon <- sf::st_sf(.all = "all", geometry = sf::st_union(sf::st_geometry(sum_polygon)))
  } else {
    grp <- sum_polygon_groups
    sum_polygon <- sum_polygon[!is.na(sum_polygon[[grp]]), ]
    if (length(class_vals) > 0) sum_polygon <- sum_polygon[sum_polygon[[grp]] %in% class_vals, ]
    sum_polygon <- sum_polygon %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
      dplyr::summarise(.groups = "drop") # summarise() on sf unions the geometries of each group
  }
  all_groups <- sort(unique(sum_polygon[[grp]]))
  ids <- conservation_areas_sf[[conservation_areas_id]]

  # 2. Intersect all conservation areas with all groups at once
  ca <- conservation_areas_sf[conservation_areas_id]
  sf::st_agr(ca) <- "constant"
  sf::st_agr(sum_polygon) <- "constant"
  inter <- sf::st_intersection(ca, sum_polygon)

  # 3. Area in km2 (assumes units of m). Touching polygons give zero area and are dropped.
  out_df <- dplyr::tibble(
      id = inter[[conservation_areas_id]],
      group = inter[[grp]],
      area_km2 = as.numeric(sf::st_area(inter)) / 1000000) %>%
    dplyr::group_by(.data$id, .data$group) %>%
    dplyr::summarise(area_km2 = round(sum(.data$area_km2), 3), .groups = "drop") %>%
    dplyr::filter(.data$area_km2 > 0)

  # 4. Optionally add a zero for every conservation area x group with no overlap
  if (fill_zeros) {
    out_df <- tidyr::expand_grid(id = ids, group = all_groups) %>%
      dplyr::left_join(out_df, by = c("id", "group")) %>%
      dplyr::mutate(area_km2 = dplyr::coalesce(.data$area_km2, 0))
  }

  # Keep the order of conservation_areas_sf, then group
  out_df <- out_df[order(match(out_df$id, ids), out_df$group), ]

  # Restore the caller's column names; drop the group column when no groups were requested
  if (is.null(sum_polygon_groups)) {
    out_df$group <- NULL
    names(out_df)[names(out_df) == "id"] <- conservation_areas_id
  } else {
    names(out_df)[names(out_df) == "id"] <- conservation_areas_id
    names(out_df)[names(out_df) == "group"] <- sum_polygon_groups
  }
  out_df
}