### seeds ###
#
#' Create a seeds table.
#'
#' The seeds table is an input to BUILDER that lists seed catchments and area targets. Seed catchments are the
#' starting catchments for constructing conservation areas. This function provides  methods for filtering the
#' catchments dataset to create seeds. Each seed requires an area target detailing the minimum area in m2 for
#' the conservation area built from that seed. Various methods are provided for adding an area target for each seed catchment.
#'
#' ### Filtering
#' If no filter arguments are provided, all catchments in \code{catchments_sf} will be added to the seeds table.
#'  \itemize{
#'    \item{To filter based on a column in \code{catchments_sf}, we recommend using dplyr to filter
#'    \code{catchments_sf} before passing it to \code{seeds()}.}
#'    \item{A \code{filter_polygon} can also be provided to filter the catchments based on their centroid (more
#'    specifically, \code{sf::st_point_on_surface()}) falling inside the \code{filter_polygon}.}
#'    }
#' ### Area targets
#' Area targets should be provided in m2.
#' Arguments to satisfy one of the following area target methods must be provided. If multiple arguments are provided, only one
#' will be used, with the priority following the order below:
#' \itemize{
#'   \item{Single-value. If \code{areatarget_value} is provided, all seeds will be assigned the area target value.}
#'   \item{Column. If \code{areatarget_col} is provided, area targets will be pulled from the column of that
#'   name in \code{catcments_sf}.}
#'   \item{Spatial join. An \code{areatarget_polygon} will be spatially joined to the catchments. Area targets
#'    will be extracted from the \code{areatarget_polygon_col} column in \code{areatarget_polygon}. If multiple polygons in
#'    \code{areatarget_polygon} intersect a catchment, the value in \code{areatarget_polygon_col} with the most overlap will be
#'    selected.}
#'   }
#'
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM.
#' @param filter_polygon Optional sf polygon object to filter catchments. Only catchments inside the polygon are kept.
#' @param areatarget_value Optional single area target value to apply to all seeds.
#' @param areatarget_col Optional column in \code{catchments_sf} holding area target values.
#' @param areatarget_polygon Optional sf polygon object holding area target values that are applied to seeds using a spatial join.
#' @param areatarget_polygon_col If \code{areatarget_polygon} provided, the column containing the area target values to join to the catchments.
#'
#' @return A tibble of seed catchments and their area targets.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Use all catchments as seeds, single area target
#' seeds(catchments_sf = builder_catchments_sample, areatarget_value = 1000000000)
#'
#' # Use all catchments as seeds, column area target
#' builder_catchments_sample$area_target <- 1000000000
#' seeds(catchments_sf = builder_catchments_sample, areatarget_col = "area_target")
#'
#'
#' # Filter based on intactness
#' library(dplyr)
#' builder_catchments_sample %>%
#'       filter(intact == 1) %>%
#'       seeds(catchments_sf = ., areatarget_value = 1000000000)
#'
#'library(sf)
#' # Filter based on polygon, assign area target via spatial join
#' ref_poly <- data.frame(
#'              lon = c(-138.4, -138.1, -138.1, -138.4,  -138.1, -138.1, -138, -138),
#'              lat = c(64.3, 64.3, 64.1, 64.1,  64.3, 64.1, 64.1, 64.3),
#'              Areatarget = c(rep(1000000000, 4),rep(2000000000,4))) %>%
#'            st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#'            group_by(Areatarget) %>%
#'            summarise(geometry = st_combine(geometry)) %>%
#'            st_cast("POLYGON") %>%
#'            st_transform(st_crs(builder_catchments_sample))
#' seeds(catchments_sf = builder_catchments_sample,
#'       filter_polygon = ref_poly,
#'       areatarget_polygon = ref_poly, areatarget_polygon_col = "Areatarget")
#'
seeds <- function(catchments_sf, filter_polygon = NULL, areatarget_value = NULL, areatarget_col = NULL, areatarget_polygon = NULL, areatarget_polygon_col = NULL){

  # SET UP
  # determine area target method.
  # priority: single value > column > polygon
  areatarget_method <- ""
  if(!is.null(areatarget_value) & is.numeric(areatarget_value)){

    areatarget_method <- "single_value"

  } else if(!is.null(areatarget_col)){
    if(areatarget_col %in% colnames(catchments_sf)){

      areatarget_method <- "column"

    } else{
      warning(paste0("areatarget_col provided but '", areatarget_col, "' not in catchments_sf"))
    }
  } else if(!is.null(areatarget_polygon)){
    if(areatarget_polygon_col %in% colnames(areatarget_polygon)){

      areatarget_method <- "polygon"

    } else{
      warning(paste0("areatarget_polygon provided but '", areatarget_polygon_col, "' not in areatarget_polygon"))
    }
  }
  if(areatarget_method == ""){
    stop("No valid area target method provided")
  }

  # CHECKS
  # if filter_polygon is provided, check for geometry in catchments and polygon
  if(!is.null(filter_polygon)){
    check_for_geometry(filter_polygon)
    check_for_geometry(catchments_sf)
  }

  # if area target method is polygon, check for geometry in catchments and polygon
  if(areatarget_method == "polygon"){
    check_for_geometry(catchments_sf)
    check_for_geometry(areatarget_polygon)
  }

  # check catchnum and convert to integer
  catchments_sf <- make_all_integer(catchments_sf, cols = "CATCHNUM")

  # FILTER
  filtered_catchments <- catchments_sf
  sf::st_agr(filtered_catchments) = "constant"

  # filter by polygon
  if(!is.null(filter_polygon)){

    catchnum_indexes <- filtered_catchments %>%
      sf::st_point_on_surface() %>% # get catchment centroids.
      sf::st_within(filter_polygon) %>% # test within for all centroids, returns a row for each match. row.id is a catchnum index, col.id is a PA index.
      as.data.frame() %>%
      dplyr::pull(row.id)

    filtered_catchments <- filtered_catchments[catchnum_indexes,] # filter using indexes from st_within
  }

  # Error if no catchments selected
  if(nrow(filtered_catchments) == 0){
    stop("No catchments selected")
  }

  # AREA TARGET
  if(areatarget_method == "single_value"){

    filtered_catchments$Areatarget <- as.integer(areatarget_value)

  } else if(areatarget_method == "column"){

    filtered_catchments$Areatarget <- as.integer(filtered_catchments[[areatarget_col]])

  } else if(areatarget_method == "polygon"){

    # Check all area targets are valid numerics
    if(!all(!is.na(suppressWarnings(as.numeric(areatarget_polygon[[areatarget_polygon_col]]))))){
      stop(paste0("All '", areatarget_polygon_col, "' values in areatarget_polygon must be numeric"))
    }

    # subset areatarget_polygon to drop any columns that might interfere with join (e.g. a CATCHNUM col)
    areatarget_polygon <- areatarget_polygon %>%
      dplyr::select(any_of(areatarget_polygon_col))

    sf::st_agr(areatarget_polygon) = "constant"
    sf::st_agr(filtered_catchments) = "constant"
    
    filtered_catchments <- suppressWarnings(filtered_catchments %>%
      sf::st_join(areatarget_polygon, left = FALSE, largest = TRUE) %>% # inner join, assign catchnum area target to polygon value with largest overlap
      dplyr::mutate(Areatarget = as.integer(ceiling(.[[areatarget_polygon_col]]))) # round up to next m2
    )

  } else{
    stop("Catchments filtered but no valid area target method provided")
  }

  out_tab <- filtered_catchments %>%
    sf::st_drop_geometry() %>%
    dplyr::as_tibble() %>%
    dplyr::select(CATCHNUM, Areatarget)

  return(out_tab)
}
