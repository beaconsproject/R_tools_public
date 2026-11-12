### catchnums_in_polygon ###
#'
#' Convert catchments to centroids and intersects with PAs
#' Output is a df of catchments in each polygon, where column names are unique ids from conservation_areas_sf. 
#' Long table would be better but this matches the output from BUILDER.

#' Using centroids can result in the inclusion of catchments that do not actually intersect conservation_areas_sf.
#' Instead use st_PointOnSurface, the calculation for which is explained here: https://gis.stackexchange.com/questions/76498/how-is-st-pointonsurface-calculated?newreg=839f6b14ce6b40c9b63956954a3ab969
#' This is essential when PAs were constructed using catchments, should be an acceptable method when non-catchment PAs are used too.
catchnums_in_polygon <- function(conservation_areas_sf, conservation_areas_id, catchments_sf){
  
  if(conservation_areas_id == "CATCHNUM"){
    conservation_areas_id <- "unique_id"
    names(conservation_areas_sf)[names(conservation_areas_sf) == "CATCHNUM"] <- "unique_id"
  }
  
  sf::st_agr(catchments_sf) = "constant"
  
  catch_within <- catchments_sf %>% 
    sf::st_point_on_surface() %>% # get catchment centroids. 
    sf::st_within(conservation_areas_sf) %>% # test within for all centroids in all conservation areas, returns a row for each match. row.id is a catchnum index, col.id is a PA index.
    as.data.frame()
  
  catchments_sf$key <- 1:nrow(catchments_sf) # add a key column to sf table. Must be an index to match st_within output
  sf_catch_key <- sf::st_drop_geometry(catchments_sf[c("key","CATCHNUM")])
  conservation_areas_sf$key <- 1:nrow(conservation_areas_sf) # add key to pa_sf
  CAs_key <- sf::st_drop_geometry(conservation_areas_sf[c("key",conservation_areas_id)])
  
  # convert indexes from st_within to catchnums using the keys to join
  tbl_long <- catch_within %>%
    dplyr::left_join(sf_catch_key, by = c("row.id" = "key")) %>%
    dplyr::left_join(CAs_key, by = c("col.id" = "key")) %>%
    dplyr::select(all_of(c("CATCHNUM", conservation_areas_id))) %>%
    dplyr::arrange(.data[[conservation_areas_id]])
  
  # convert long table to wide table with missing values as NA
  out_tab <- long_to_wide(tbl_long, conservation_areas_id, "CATCHNUM")
  return(out_tab)
}
