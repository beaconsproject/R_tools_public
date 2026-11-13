### build_network_polygons ###
#
#' Combine a set of conservation area polygons into a network.
#'
#' Takes a set of individual conservation area polygons and/or generic reserve polygons and combines them into network 
#' polygons based on a list of network names. Generally used for combining conservation areas from the beaconsbuilder package
#' into networks, but can also be used to build networks using any input reserve polygons (e.g. existing protected areas).
#' 
#' For large lists of networks (>10,000), we recommend subsetting the network_list and making multiple calls
#' to \code{\link{build_network_polygons}}.
#'
#' @param conservation_areas_sf sf object with unique id column named \code{network}, typically the output from
#'  [dissolve_catchments_from_table()].
#' @param conservation_areas_id character specify unique identifier. Default is network.
#' @param network_list Vector of strings detailing network names to be built. Typically network names built using [gen_network_names()] 
#' (e.g. PB_0001__PB_0002). All network names must include names from the conservation_areas_sf \code{network} column, separated by \code{"__"} (double underscore).
#'
#' @return A sf object of networks.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample, 
#'   "network",
#'   dissolve_list = colnames(builder_table_sample))
#' network_names <- gen_network_names(conservation_areas$network, 2)
#' build_network_polygons(conservation_areas, network_names)
#' 
#' build_network_polygons(conservation_areas, c("PB_0001__PB_0002"))

build_network_polygons <- function(conservation_areas_sf, conservation_areas_id = "network", network_list){
  
  missing_cols <- check_colnames(conservation_areas_sf, conservation_areas_id)
  if (length(missing_cols) > 0) stop("Column(s) missing from conservation_areas_sf: ", paste(missing_cols, collapse = ", "))
  check_for_geometry(conservation_areas_sf)
  
  my_list <- list()
  for(net in network_list){
    
    # get benchmark names and check they are in conservation_areas_sf$network
    nets <- sep_network_names(net)
    check_conservation_area_names(nets, conservation_areas_sf)
    
    # dissolve the network
    row_sfc <- conservation_areas_sf %>%
      dplyr::filter(.data$network %in% nets) %>%
      sf::st_union()
    row <- sf::st_sf(network = net, geometry = row_sfc)
    
    my_list <- append(my_list, list(row))
  }
  
  if(utils::packageVersion("dplyr") > "1.0.0"){
    out_sf <- dplyr::bind_rows(my_list) # this is much faster than do.call but requires dplyr >0.9 for bind_rows to work with sf objects
  } else{
    out_sf <- do.call(rbind, my_list)
  }
  return(out_sf)
}
