### append_reserve ###
#
#' Adds a polygon to a sf object.
#'
#' Used to add reserves that are not output from the beaconsbuilder package to a sf object of conservation areas.
#' For example, adding an existing protected area (PA) to an sf object of conservation_areas build using the 
#' beaconsbuilder package. This allows the appended reserve to be included in networks built using [build_network_polygons()].
#'
#' @param conservation_areas_sf sf object with unique id column named \code{network}, typically the output from
#'  [dissolve_catchments_from_table()].
#'@param conservation_areas_id
#' @param add_reserve sf object to add as a single additional reserve to conservation_areas_sf. All features will
#' be dissolved into a single POLYGON or MULTIPOLYGON feature to append to conservation_areas_sf. If multiple add_reserve
#' object are required, add them using multiple calls to [append_reserve()].
#' @param reserve_name String that will become the reserve name in the conservation_areas_sf \code{network} column.
#'
#' @return sf object matching conservation_areas_sf.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library(sf)
#' library(dplyr)
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample, 
#'   "network",
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' pa_1 <- data.frame(lon = c(-85, -82.5, -83), lat = c(51, 51.5, 50.5)) %>%
#'           sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#'           dplyr::summarise(geometry = st_combine(geometry)) %>%
#'           sf::st_cast("POLYGON") %>%
#'           sf::st_transform(st_crs(conservation_areas)) %>%
#'           dplyr::mutate(pa_name="a protected area", ha=9999, class="provincial")
#' append_reserve(conservation_areas, pa_1, "PA_1")
append_reserve <- function(conservation_areas_sf, conservation_areas_id, add_reserve, reserve_name){
  
  conservation_areas_sf <- check_colnames(conservation_areas_sf, cols = conservation_areas_id)
  check_for_geometry(conservation_areas_sf)
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(add_reserve))
  
  # dissolve add_reserve into one feature
  add_reserve <- add_reserve %>%
    dplyr::summarise(geometry = sf::st_union(.data$geometry)) %>%
    dplyr::mutate(network = reserve_name)
  
  # append to conservation_areas_sf
  conservation_areas_sf <- rbind(conservation_areas_sf, add_reserve)
  
  return(conservation_areas_sf)
}
