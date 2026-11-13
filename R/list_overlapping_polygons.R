### list_overlapping_polygons ###
#
#' List all pairs of touching/overlapping conservation areas.
#'
#' For all touching/overlapping polygons in the input file, returns the network_name of the overlapping
#' pairs.
#' 
#' This can be converted to a list of pairs of individual names using \code{sep_network_names(list_overlapping_polygons())}. 
#' 
#'
#' @param conservation_areas_sf sf object with unique id column named \code{network}, typically the output from
#'  [dissolve_catchments_from_table()], possibly with additional reserves added using [append_reserve()].
#' @param conservation_areas_id character specify unique identifier. Default is network.
#' @return Vector of network names constructed using the \code{"__"} separator.
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample, 
#'   "network",
#'   dissolve_list = colnames(builder_table_sample))
#' list_overlapping_polygons(conservation_areas)
#' 
#' sep_network_names(list_overlapping_polygons(conservation_areas))

list_overlapping_polygons <- function(conservation_areas_sf, conservation_areas_id = "network"){
  
  missing_cols <- check_colnames(conservation_areas_sf, conservation_areas_id)
  if (length(missing_cols) > 0) stop("Column(s) missing from conservation_areas_sf: ", paste(missing_cols, collapse = ", "))
  
  df <- as.data.frame(sf::st_intersects(conservation_areas_sf, conservation_areas_sf)) # get pairwise intersects
  
  # convert to list, drop duplicates and self intersects, sort names
  l <- unique(
    lapply(1:nrow(df), function(x){ 
      if(df[x,1] != df[x,2]){
        sort(c(df[x,1], df[x,2]))
        }
      })
    )
  l <- l[!sapply(l, is.null)] # remove null
  
  return_names <- sapply(l, function(x){ # convert indexes to names
    paste0(sort(c(conservation_areas_sf$network[x[1]], conservation_areas_sf$network[x[2]])), collapse = "__")
  })
  return(return_names)
}
