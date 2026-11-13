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
#'        Options are `"CENTROID"` (default, catchments whose centroid falls inside the conservation area) or
#'        `"INTERSECT"` (catchments touching the conservation area).
#'
#' @return A tibble of reserve seed catchments and their associated area targets.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' 
#' # Conservation areas built from the sample builder table
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, builder_table_sample, "network",
#'   dissolve_list = c("PB_0001", "PB_0002"))
#'
#' reserve_seeds(catchments_sf = catchments_sample,
#'       conservation_areas_sf = conservation_areas, conservation_areas_id = "network",
#'       areatarget_value = 1000000000, joinType = "CENTROID")
#'
#' # Also include catchments touching the edge of the conservation areas
#' reserve_seeds(catchments_sample, conservation_areas, "network",
#'       areatarget_value = 1000000000, joinType = "INTERSECT")
#'
reserve_seeds <- function (catchments_sf, conservation_areas_sf, conservation_areas_id, areatarget_value= NULL, joinType = c("CENTROID", "INTERSECT")){
  joinType <- match.arg(joinType)

  # avoid scientific notation when pasting CATCHNUMs, restore user setting on exit
  old_opts <- options(scipen = 999)
  on.exit(options(old_opts), add = TRUE)

  # JOIN TYPE
  if(joinType =="INTERSECT"){
    catch_join <- suppressWarnings(st_intersects(conservation_areas_sf, catchments_sf))
  }else if(joinType =="CENTROID"){
    catch_centroids <- suppressWarnings(st_centroid(catchments_sf))
    catch_join <- suppressWarnings(st_intersects(conservation_areas_sf, catch_centroids, sparse = TRUE))
  }
  
  out_tab <- tibble::tibble(
    !!conservation_areas_id := conservation_areas_sf[[conservation_areas_id]],
    Areatarget = if (!is.null(areatarget_value)) as.integer(areatarget_value), # integer, same as seeds(); NULL drops the column
    CATCHNUM = sapply(catch_join, function(idx) paste(catchments_sf$CATCHNUM[idx], collapse = ","))
  )
  
  return(out_tab)
}
