### reserve_seeds ###
#
#' Create a seeds table based on conservation areas.
#'
#' This function identifies as seeds catchments that intersect conservation areas.
#' Each seed requires an area target detailing the minimum area in m2 for the conservation area built from that seed. 
#' 
#' @param catchments_sf An `sf` object of the catchments dataset. Must contain a unique identifier column named `CATCHNUM`.
#' @param conservation_areas_sf An `sf` polygon object representing conservation areas from which the catchments will be extracted.
#' @param conservation_areas_id A character string. Name of the column holding conservation area unique identifiers. 
#' @param areatarget_value Optional. A single numeric value specifying the area target to apply to all seeds.
#' @param joinType A character string. Type of spatial join to identify catchments within conservation areas. 
#'        Options are `"INTERSECT"` or `"CENTROID"`. Default is `"CENTROID"`.
#'
#' @return A tibble of reserve seed catchments and their associated area targets.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' 
#'library(sf)
#' # Filter based on polygon, assign area target via spatial join
#' # Create a polygon from coordinates
#'CAs_coords <- data.frame(
#'  lon = c(-2014941, -2014941, -2004941, -2004941, -2014941),  # repeat first point to close the polygon
#'  lat = c(802374.6, 812374.6, 812374.6, 802374.6, 802374.6)
#')
#'
#'# Turn into an sf polygon
#'ref_poly <- st_sf(
#'  ID = "CA_01",
#'  geometry = st_sfc(
#'    st_polygon(list(as.matrix(CAs_coords)))
#'  ),
#'  crs = st_crs(builder_catchments_sample)  # use same CRS as your other data
#') 
#' reserve_seeds(catchments_sf = builder_catchments_sample,
#'       conservation_areas_sf = ref_poly, conservation_areas_id = "ID",
#'       areatarget_value = 10000, joinType = "CENTROID")
#'
reserve_seeds <- function (catchments_sf, conservation_areas_sf, conservation_areas_id, areatarget_value= NULL, joinType = "CENTROID"){
  options(scipen = 999)
  
  # JOIN TYPE
  if(joinType =="INTERSECT"){
    catch_join <- suppressWarnings(st_intersects(conservation_areas_sf, catchments_sf))
  }else if(joinType =="CENTROID"){
    catch_centroids <- suppressWarnings(st_centroid(catchments_sf))
    catch_join <- suppressWarnings(st_intersects(conservation_areas_sf, catch_centroids, sparse = TRUE))
  }
  
  out_tab <- tibble::tibble(
    !!conservation_areas_id := conservation_areas_sf[[conservation_areas_id]],
    Areatarget = areatarget_value,
    CATCHNUM = sapply(catch_join, function(idx) paste(catchments_sf$CATCHNUM[idx], collapse = ","))
  )
  
  return(out_tab)
}
