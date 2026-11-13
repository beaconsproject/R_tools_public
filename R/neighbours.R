### neighbours ###
#
#' Create a neighbours table listing neighbours for each catchment.
#'
#' For an sf object of catchments with unique CATCHNUM id's, calculates a list of neighbouring CATCHNUM pairs
#' and returns them in a long tibble. Neighbours are defined as having at least on point in common (within 0.1m).
#'
#'
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM, in a projected CRS (metres).
#'
#' @return A tibble of neighbouring pairs with columns \code{CATCHNUM} and \code{neighbours}.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' neighbours(catchments_sample)
#'
neighbours <- function(catchments_sf){
  # Generate neighbours where single point is shared
  # this was originally done in Python using the GenerateSpatialWeightsMatrix() function
  # here we use the SF package
  # pattern = "****T****" in the st_relate() function matches any intersecting polygons.
  # More info on st_relate at:
  # https://www.rdocumentation.org/packages/sf/versions/0.7-7/topics/st_relate
  # Queen pattern found here: https://github.com/r-spatial/sf/issues/234

  # the 0.1 buffer below is in map units, so the CRS must be projected (metres), not lat/long degrees
  if(isTRUE(sf::st_is_longlat(catchments_sf))){
    stop("catchments_sf must use a projected CRS in metres, not longitude/latitude")
  }

  # check catchnum and convert to integer
  catchments_sf <- make_all_integer(catchments_sf, cols= "CATCHNUM")

  st_queen <- function(a, b = a) sf::st_relate(a, b, pattern = "****T****") # this tests for an intersect of at least one shared point between a and b
  nbr_df <- as.data.frame(st_queen(sf::st_buffer(catchments_sf, dist=0.1)))

  # replace index values with catchnum values using a key
  catchments_sf$key <- 1:nrow(catchments_sf) # add a key column to sf table. Must be an index so it matches the index assigned to the NB_QUEEN column
  sf_catch_key <- sf::st_drop_geometry(catchments_sf[c("key","CATCHNUM")])

  nbr_df <- nbr_df %>%
    dplyr::left_join(sf_catch_key, by = c("row.id" = "key")) %>%
    dplyr::left_join(sf_catch_key, by = c("col.id" = "key")) %>%
    dplyr::select(CATCHNUM.x, CATCHNUM.y)

  names(nbr_df) <- c("CATCHNUM", "neighbours") # rename output columns

  # remove cases where CATCHNUM is its own NEIGHBOUR
  nbr_df <- nbr_df %>%
    dplyr::filter(CATCHNUM != neighbours) %>%
    dplyr::as_tibble()

  # add key needed by BUILDER
  nbr_df$key <- as.integer(0:(nrow(nbr_df)-1))

  return(nbr_df)
}
