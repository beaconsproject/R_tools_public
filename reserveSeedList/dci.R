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
#' @param conservation_area_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param buffer_width Width of buffer to apply to stream segments. Defaults to 0.1. Used to
#' ensure adjacent stream segments are connected during analysis.
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
calc_dci <- function(conservation_area_sf, stream_sf, column_id = network, buffer_width = 0.1){
  #browser()
  stopifnot(sf::st_crs(conservation_area_sf) == sf::st_crs(stream_sf))

  # set stream and conservation area attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(conservation_area_sf) = "constant"
  
  # clip streams to full area of all conservation areas
  reserve_dci <- conservation_area_sf %>%
    sf::st_intersection(stream_sf, sf::st_union()) %>% # get streams just for the required region
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::summarise(geometry = sf::st_union(geometry)) %>% # merge into single feature
    sf::st_intersection(conservation_area_sf) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select({{column_id}}) %>% # drop attributes except network id
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(column_id)) %>% # for each network...
    dplyr::summarise(L = sum(stream_length), dci = sum((stream_length*stream_length) / (L*L)))  # calculate L2 then use to calculate dci
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  dci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_dci, by = column_id) %>%
    tidyr::replace_na(list(dci=0)) %>%
    dplyr::pull(dci) %>%
    round(3)
  
  return(dci)
}

#' calc_lwdci
#' 
#' Calculate stream network length-weighted dendritic connectivity index (DCI) for a set of input polygons.
#'
#' Dendritic connectivity measures the longitudinal connectivity of a river network inside a set of conservation areas
#' based on the expected probability of an organism being able to move between two random points in the
#' network [(Cote et al 2009)](https://link.springer.com/article/10.1007/s10980-008-9283-y).
#'
#' DCI = sum ( li2 / L2 )
#' li - length of stream section
#' L - total length of all stream sections in conservation area
#' 
#' Values range between 0 and 1, with 1 indicating a completely connected river network.
#' 
#' For length weighted DCI, a separate DCI measure is calculated for each group of BASIN streams 
#' within the conservation area, then a weighted average is calculated where the weights are the lengths 
#' of streams in each BASIN.
#'
#' @param conservation_area_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param buffer_width Width of buffer to apply to stream segments. Defaults to 0.1. Used to
#' ensure adjacent stream segments are connected during analysis.
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
#' calc_lwdci(conservation_areas, streams_sample)
#' 

calc_lwdci <- function(conservation_area_sf, stream_sf, column_id = network, buffer_width = 0.1){
  
  stopifnot(sf::st_crs(conservation_area_sf) == sf::st_crs(stream_sf))
  
  # set stream and reserve attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(conservation_area_sf) = "constant"
  
  # clip streams to full area of all reserves
  stream_prepped <- 
    sf::st_intersection(stream_sf, sf::st_union(conservation_area_sf$geometry)) %>%
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::group_by(BASIN) %>%
    dplyr::summarise(geometry = sf::st_union(geometry)) # merge into single feature per BASIN
  
  sf::st_agr(stream_prepped) = "constant"
  
  reserve_lwdci <- sf::st_intersection(conservation_area_sf, stream_prepped) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select({{column_id}}, BASIN) %>% # drop attributes except network and BASIN
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(column_id), BASIN) %>% # for each network and basin...
    dplyr::summarise(L_basin = sum(stream_length), 
                     dci_basin = sum((stream_length*stream_length) / (L_basin*L_basin)), .groups = "drop") %>% # calculate L2 then use to calculate dci
    dplyr::group_by(!!sym(column_id)) %>%
    dplyr::summarise(lwdci = sum((L_basin / sum(L_basin)) * dci_basin), .groups = "drop") # calc lwdci. Weight is proportion of reserve streams in each basin
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  lwdci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_lwdci, by = column_id) %>%
    tidyr::replace_na(list(lwdci=0)) %>%
    dplyr::pull(lwdci) %>%
    round(3)
  
  return(lwdci)
}
