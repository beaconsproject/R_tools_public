### calc_dci ###
#
#' Calculate dendritic connectivity index (DCI) for a set of input polygons.
#' 
#' Values range between 0 and 1, with 1 indicating a completely connected river network.
#'
#' Dendritic connectivity measures the longitudinal connectivity of a river network inside a set of reserves
#' based on the expected probability of an organism being able to move between two random points in the
#' network [(Cote et al 2009)](https://link.springer.com/article/10.1007/s10980-008-9283-y).
#'
#' DCI = sum ( li2 / L2 )
#' li - length of stream section
#' L - total length of all stream sections in conservation area
#' 
#'
#' @param conservation_areas_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param conservation_areas_id character specify unique identifier. Default is network.
#' @param buffer_width Width of buffer to apply to stream segments. Defaults to 0.1. Used to
#' ensure adjacent stream segments are connected during analysis.
#' @param progress Optional callback function for reporting progress. The function
#'   should accept two arguments: `value`, a numeric value between 0 and 1
#'   indicating the progress, and `detail`, a character string describing the
#'   current processing step. Defaults to `NULL`.
#'
#' @return Vector of numeric DCI values matching the input features. 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' calc_dci(conservation_areas, streams_sample)
#' 
calc_dci <- function(conservation_areas_sf, stream_sf, conservation_areas_id = "network", buffer_width = 0.1, progress = NULL){
  
  report <- function(value, detail) {
    if (!is.null(progress)) {
      progress(value, detail)
    }
  }
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(stream_sf))
  
  # set stream and conservation area attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(conservation_areas_sf) = "constant"
  
  report(0.05, "Clipping streams to conservation areas")
  
  conservation_area_union <- sf::st_union(sf::st_geometry(conservation_areas_sf))
  
  stream_candidates <- stream_sf |>
    dplyr::filter(BASIN != -1) |>
    sf::st_filter(conservation_area_union)
  
  stream_prepped <- sf::st_intersection(
    stream_candidates,
    conservation_area_union
  )
  
  report(0.25, "Buffering stream sections")
  buffered_streams <- stream_prepped |>
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND")
  
  report(0.50, "Merging connected stream sections")
  merged_streams <- sf::st_sf(
     geometry = sf::st_union(sf::st_geometry(buffered_streams))
  )
  
  report(0.70, "Calculating components within each conservation area")
  reserve_dci <- sf::st_intersection(merged_streams, conservation_areas_sf) |>
    dplyr::select(dplyr::all_of(conservation_areas_id)) |>
    sf::st_cast("MULTIPOLYGON", warn = FALSE) |>
    sf::st_cast("POLYGON", warn = FALSE) |>
    dplyr::mutate(
      stream_length = as.numeric(sf::st_area(geometry)) / buffer_width
    ) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(!!sym(conservation_areas_id)) |>
    dplyr::summarise(
      dci = sum(stream_length^2) / sum(stream_length)^2,
      .groups = "drop"
    )
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  report(0.90, "Joining DCI values to conservation areas")
  dci <- dplyr::left_join(sf::st_drop_geometry(conservation_areas_sf), reserve_dci, by = conservation_areas_id) %>%
    dplyr::mutate(dci = tidyr::replace_na(dci, 0)) %>%
    dplyr::pull(dci) %>%
    round(3)
  
  report(1, "Complete")
  return(dci)
}