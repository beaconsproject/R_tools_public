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
  catchments_sf <- make_catchnum_integer(catchments_sf)

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
      dplyr::select(areatarget_polygon_col)

    sf::st_agr(areatarget_polygon) = "constant"
    sf::st_agr(filtered_catchments) = "constant"
    browser()
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
neighbours <- function(catchments_sf){
  # Generate neighbours where single point is shared
  # this was originally done in Python using the GenerateSpatialWeightsMatrix() function
  # here we use the SF package
  # pattern = "****T****" in the st_relate() function matches any intersecting polygons.
  # More info on st_relate at:
  # https://www.rdocumentation.org/packages/sf/versions/0.7-7/topics/st_relate
  # Queen pattern found here: https://github.com/r-spatial/sf/issues/234

  # check catchnum and convert to integer
  catchments_sf <- make_catchnum_integer(catchments_sf)

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
make_catchnum_integer <- function(catchments_sf){

  # check CATCHNUM exists
  if(!"CATCHNUM" %in% names(catchments_sf)){
    stop("Catchments must contain column 'CATCHNUM'")
  }

  catchments_sf$CATCHNUM <- as.integer(catchments_sf$CATCHNUM)

  return(catchments_sf)
}
check_for_geometry <- function(in_sf){

  if(!"geometry" %in% names(in_sf)){
    stop("Must contain column: geometry")
  }
}
