### neighbours ###
#
#' Create a neighbours table listing neighbours for each catchment.
#'
#' For an sf object of catchments with unique CATCHNUM id's, calculates a list of neighbouring CATCHNUM pairs
#' and returns them in a long tibble. Neighbours are defined as having at least on point in common (within 0.1m).
#'
#'
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM .
#'
#' @return A tibble of neighbouring pairs with columns \code{CATCHNUM} and \code{neighbours}.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' neighbours(builder_catchments_sample)
#'
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

### reserve_seeds ###
#
#' Create a seeds table based on conservation areas.
#'
#' This function identifies as seeds catchments that intersect conservation areas.
#' Each seed requires an area target detailing the minimum area in m2 for the conservation area built from that seed. 
#' 
#' @param catchments_sf An `sf` object of the catchments dataset. Must contain a unique identifier column named `CATCHNUM`.
#' @param CAs_sf An `sf` polygon object representing conservation areas from which the catchments will be extracted.
#' @param CAs_name A character string. Name of the column holding conservation area unique identifiers. 
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
#'       CAs_sf = ref_poly, CAs_name = "ID",
#'       areatarget_value = 10000, joinType = "CENTROID")
#'
reserve_seeds <- function (catchments_sf, CAs_sf, CAs_name, areatarget_value= NULL, joinType = "CENTROID"){
  options(scipen = 999)
  
  # JOIN TYPE
  if(joinType =="INTERSECT"){
    catch_join <- suppressWarnings(st_intersects(CAs_sf, catchments_sf))
  }else if(joinType =="CENTROID"){
    catch_centroids <- suppressWarnings(st_centroid(catchments_sf))
    catch_join <- suppressWarnings(st_intersects(CAs_sf, catch_centroids, sparse = TRUE))
  }
  
  out_tab <- tibble::tibble(
    !!CAs_name := CAs_sf[[CAs_name]],
    Areatarget = areatarget_value,
    CATCHNUM = sapply(catch_join, function(idx) paste(catchments_sf$CATCHNUM[idx], collapse = ","))
  )
  
  return(out_tab)
}

### builder ###
#
#' Call the BUILDER software from R.
#'
#' Prepares and passes all input tables and variables to the BUIDLER executable and returns the builder table to R.
#'
#' This function prepares the input tables (neighbours, seeds, catchments) required by BUILDER in a temp file, then creates
#' a string of parameters that is passed to the BUILDER executable using \code{system()}. The BUILDER output tables are saved to the
#' temp file, and the output table listing conservation area names and catchment lists is returned to the R session.
#'
#' This function only returns the BUILDER output table describing the lists of catchments making up each conservation area. We typically recommend
#' using the \code{beaconstools} package to calculate hydrology metrics, but hydrology information as well as other summary information is
#' produced by BUILDER.
#' If users want access to the complete set of BUILDER output tables (e.g. hydrology and summary information), a valid output directory
#' can be provided in which case all input and output tables will be saved to that directory. If no output directory is provided a temp
#' folder will be used and deleted after the function completes. As well as providing access to the full list of BUILDER output tables,
#' an output directory is useful in a looped workflow with multiple calls to BUILDER producing multiple sets of output files. In this
#' workflow a different output directory can be used to save each set of output files which can later be processed using the \code{beaconstools} package.
#' The main output file listing conservation area names and their catchments is named "*_COLUMN_All_Unique_BAs.csv".
#'
#'
#' The main parameters to adjust when building conservation areas are:
#' \itemize{
#'   \item{Area target: Set in the seeds table, this defines the minimum size for each conservation area.}
#'   \item{Catchment-level intactness: Only catchments with intactness values greater than this value will be added to the conservation area.}
#'   \item{Conservation area intactness: Only conservation areas with an area-weighted intactness greater than or equal to this value will be returned.}
#'   \item{Area-target proportion: The proportion of the Area target that conservation areas need to meet. Allows output of conservation areas smaller than
#'   the area target (e.g. if no conservation areas meet 100% of the area target, area target multiplier can be lowered to 0.9 to search for smaller
#'   conservation areas).}
#'   }
#'
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM.
#' @param data_source Specify the seed type (catchment or reserve). Default is catchment. 
#' @param seeds Seeds table from \code{seeds()} listing seed catchments and area targets.
#' @param reserve_name Reserve column listing the name of the reserve. Only used if data_source Reserve is selected. 
#' @param neighbours Neighbours table from \code{neighbours()} listing all neighbouring pairs of catchments.
#' @param out_dir If provided, input (seeds, neighbours and catchments) files and output BUILDER tables will be saved to this directory. Otherwise
#' a temp directory will be used. Function will attempt to create the directory if it doesn't already exist.
#' @param builder_local_path Directory path to the benchmarkbuilder executable. 
#' @param catchment_level_intactness Minimum intactness value for catchment inclusion (between 0-1). i.e. if value of 1 is used, only 100% intact
#' catchments will be used.
#' @param conservation_area_intactness Minimum area-weighted intactness of final conservation areas. Only conservation areas meeting this value will be returned.
#' @param area_target_proportion The proportion of the Area target that conservation areas need to meet (0-1) to be returned by the function.
#' @param area_type If 'land', only terrestrial areas are counted towards the area target. If 'water', only water areas are counter. If 'landwater',
#' all areas are counted.
#' @param construct_conservation_areas Should builder build conservation areas?
#' @param area_target_multiplier Multiplier applied to the area target in the seed list that allows for adjustments to the area target without having
#' to remake the seeds table.
#' @param handle_isolated_catchments Should small isolated catchments be merged into the conservation areas?
#' @param output_upstream Should upstream catchments be calculated?
#' @param output_downstream Should downstream catchments be calculated?
#' @param output_hydrology_metrics Should hydrology metrics be calculated?
#' @param area_land Catchments column listing the area of land in each catchment.
#' @param area_water Catchments column listing the area of water in each catchment.
#' @param skeluid Catchments column listing the skeluid field.
#' @param catchnum Catchments column listing the unique catchments ID.
#' @param subzone Catchments column listing the hydrological subzone.
#' @param zone Catchments column listing the hydrological zone.
#' @param basin Catchments column listing the basin field.
#' @param order1 Catchments column listing the order 1 field.
#' @param order2 Catchments column listing the order 2 field.
#' @param order3 Catchments column listing the order 3 field.
#' @param stream_length Catchments column listing stream length.
#' @param intactness Catchments column listing the intactness value.
#' @param isolated Catchments column identifying catchments isolated from the stream network.
#' @param unique_identifier Prefix for conservation area names.
#' @param handler_summary Should a summary listing counts of available conservation areas at different intactness and area target proportions be created?
#' @param summary_intactness_props If a summary is created, provide the intactness proportions to be summarised in the format '0.7,0.8,0.9'.
#' @param summary_area_target_props If a summary is created, provide the area target proportions to be summarised in the format '0.7,0.8,0.9'.
#'
#' @return A tibble with column names representing conservation areas, and rows of catchments making up the conservation areas.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' library(dplyr)
#' nghbrs <- neighbours(builder_catchments_sample)
#' # select 10 seeds with 100% intactness for the example. Build conservation areas to 500km2.
#' seed <- builder_catchments_sample %>%
#'               filter(intact == 1) %>%
#'               seeds(catchments_sf = ., areatarget_value = 500000000)
#' seed <- seed[1:10,]
#' builder(catchments_sf = builder_catchments_sample, seeds = seed, neighbours = nghbrs)
builder <- function(catchments_sf, data_source = "catchment", seeds, reserve_name= NULL, neighbours, # input tables
                    out_dir = NULL, # output dir
                    builder_local_path,
                    catchment_level_intactness = 1, conservation_area_intactness = 1, area_target_proportion = 1, # main parameters
                    area_type = "land", construct_conservation_areas = TRUE,  area_target_multiplier = 1,
                    handle_isolated_catchments = TRUE, output_upstream = FALSE, output_downstream = FALSE, output_hydrology_metrics = FALSE,
                    area_land = "Area_Land", area_water = "Area_Water",
                    skeluid = "SKELUID", catchnum = "CATCHNUM", subzone = "FDA", zone = "MDA", basin = "BASIN", order1 = "ORDER1",
                    order2 = "ORDER2", order3 = "ORDER3", stream_length = "STRMLEN", intactness = "intact", isolated = "Isolated",
                    unique_identifier = "PB", handler_summary = FALSE, summary_intactness_props = '""', summary_area_target_props = '""'
){
  
  # builder path
  #if(exists("builder_local_path")){
  builder_cmd <- file.path(builder_local_path, "BenchmarkBuilder_cmd.exe") # This is a temporary solution to avoid publishing builder.exe on github during the review process.
  #} else{
  #  builder_cmd <- system.file("exec", "BenchmarkBuilder_cmd.exe", package = "beaconsbuilder")
  #}
  
  # set up output directory. Default is a temp folder. out_dir can be used if builder output should be saved, e.g. in a looped analysis.
  if(is.null(out_dir)){
    tmpdir <- file.path(tempdir(), "builder")
    if(!dir.exists(tmpdir)){
      dir.create(tmpdir)
    }
    outdir <- tmpdir
  } else{
    outdir <- out_dir
    # check no spaces
    if(grepl(" ", outdir)){
      stop("out_dir cannot have spaces")
    }
    if(!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE)
    }
  }
  
  # Count benchmark tables in out_dir
  benchmarks_out_pre <- length(list.files(outdir, pattern = "COLUMN_All_Unique_BAs.csv"))
  
  # save seeds table to outdir
  message("checking seeds table")
  
  if(data_source == "catchment"){
    seeds <- make_all_integer(seeds)
    check_colnames(seeds, "seeds", c("CATCHNUM", "Areatarget"))
  }else if(data_source == "Reserve"){
    check_colnames(seeds, "seeds", c(reserve_name, "Areatarget", "CATCHNUM"))
  }else{
    stop("data_source must be either catchment or Reserve")
  }
  
  check_seeds_areatargets(seeds)
  check_seeds_in_catchments(seeds, catchments_sf)
  utils::write.csv(seeds, file.path(outdir, "seeds.txt"), row.names = FALSE)
  
  # save neighbours table to outdir
  message("checking neighbours table")
  check_colnames(neighbours, "neighbours", c("CATCHNUM", "neighbours", "key"))
  neighbours <- make_all_integer(neighbours)
  utils::write.csv(neighbours, file.path(outdir, "neighbours.csv"), row.names = FALSE)
  
  # save catchment attributes
  catchments_csv <- sf::st_drop_geometry(catchments_sf)
  check_colnames(catchments_csv, "catchments_sf", c(area_land, area_water, skeluid, catchnum, subzone, zone, basin, order1, order2, order3, stream_length, intactness, isolated)) # check all provided columns are present
  catchments_csv <- make_all_character(catchments_csv, c(order1, order3, zone, subzone, basin)) # check type
  catchments_csv <- make_all_integer(catchments_csv, c(catchnum, skeluid, isolated))
  catchments_csv <- make_all_numeric(catchments_csv, c(order2, area_land, area_water, intactness))
  
  catchments_file <- file.path(outdir, "catchments.csv")
  catchments_table_name <- "catchments"
  utils::write.csv(catchments_csv, catchments_file, row.names = FALSE)
  
  # convert TRUE FALSE to 1 0
  construct_conservation_areas <- logical_to_integer(construct_conservation_areas)
  handle_isolated_catchments <- logical_to_integer(handle_isolated_catchments)
  output_upstream <- logical_to_integer(output_upstream)
  output_downstream <- logical_to_integer(output_downstream)
  output_hydrology_metrics <- logical_to_integer(output_hydrology_metrics)
  handler_summary <- logical_to_integer(handler_summary)
  
  # set system variables
  seeds_file <- file.path(outdir, "seeds.txt")
  seeds_table_name <- "seeds"
  seeds_catchnum_field <- "CATCHNUM"
  seeds_areatarget_field <- "Areatarget"
  neighbours_file <- file.path(outdir, "neighbours.csv")
  neighbours_file_name <- "neighbours"
  neighbours_catchnum_field <- "CATCHNUM"
  neighbours_neighbours_field <- "neighbours"
  
  # construct system call
  if(data_source == "catchment"){
    system_string <- paste(sep = " ",
                           builder_cmd,
                           "builder", "1",
                           outdir,
                           "catchment",
                           seeds_file,
                           seeds_table_name,
                           seeds_catchnum_field,
                           seeds_areatarget_field,
                           area_type,
                           construct_conservation_areas,
                           catchment_level_intactness,
                           area_target_multiplier,
                           handle_isolated_catchments,
                           output_upstream,
                           output_downstream,
                           output_hydrology_metrics,
                           catchments_file,
                           catchments_table_name,
                           area_land,
                           area_water,
                           skeluid,
                           catchnum,
                           subzone,
                           basin,
                           order1,
                           order2,
                           order3,
                           stream_length,
                           intactness,
                           zone,
                           isolated,
                           neighbours_file,
                           neighbours_file_name,
                           neighbours_catchnum_field,
                           neighbours_neighbours_field
    )
  }else{
    if(is.null(reserve_name)){
      stop("You must provide the column name holder the reserve name")
    }else{
      if(reserve_name %in% colnames(seeds)){
        system_string <- paste(sep = " ",
                               builder_cmd,
                               "builder", "1",
                               outdir,
                               "Reserve",
                               seeds_file,
                               seeds_table_name,
                               reserve_name,
                               seeds_catchnum_field,
                               seeds_areatarget_field,
                               area_type,
                               construct_conservation_areas,
                               catchment_level_intactness,
                               area_target_multiplier,
                               handle_isolated_catchments,
                               output_upstream,
                               output_downstream,
                               output_hydrology_metrics,
                               catchments_file,
                               catchments_table_name,
                               area_land,
                               area_water,
                               skeluid,
                               catchnum,
                               subzone,
                               basin,
                               order1,
                               order2,
                               order3,
                               stream_length,
                               intactness,
                               zone,
                               isolated,
                               neighbours_file,
                               neighbours_file_name,
                               neighbours_catchnum_field,
                               neighbours_neighbours_field
        )
      }else {
        stop("reserve_name must be contained in the reserve seed file provided")
      }
    }
  }
  
  
  # call BUILDER
  system(system_string)
  
  # prep for HANDLER
  h1 <- list.files(outdir, pattern = "_ROW.csv", full.names = TRUE)
  h1 <- ifelse(length(h1) > 0, h1[[length(h1)]], '""') # if multiple, return the last one which should be most recent
  h2 <- list.files(outdir, pattern = "_ROW_AREALAND.csv", full.names = TRUE)
  h2 <- ifelse(length(h2) > 0, h2[[length(h2)]], '""') # if multiple, return the last one which should be most recent
  h3 <- list.files(outdir, pattern = "_ROW_AREALAND.csv", full.names = TRUE)
  h3 <- ifelse(length(h3) > 0, h3[[length(h3)]], '""') # if multiple, return the last one which should be most recent
  h4 <- list.files(outdir, pattern = "_ROW_UPSTREAM_CATCHMENTS.csv", full.names = TRUE)
  h4 <- ifelse(length(h4) > 0, h4[[length(h4)]], '""') # if multiple, return the last one which should be most recent
  h5 <- list.files(outdir, pattern = "_ROW_DOWNSTREAM_CATCHMENTS.csv", full.names = TRUE)
  h5 <- ifelse(length(h5) > 0, h5[[length(h5)]], '""') # if multiple, return the last one which should be most recent
  
  # prep HANDLER system string
  handler_string <- paste(sep = " ",
                          builder_cmd,
                          "handler",
                          h1, h2, h3, h4, h5,
                          conservation_area_intactness,
                          area_target_proportion,
                          unique_identifier,
                          handler_summary,
                          summary_intactness_props,
                          summary_area_target_props
  )
  # call HANDLER
  system(handler_string)
  
  # Load output tables
  benchmarks_out <- list.files(outdir, pattern = "COLUMN_All_Unique_BAs.csv")
  if(length(benchmarks_out) > benchmarks_out_pre){ # If new table was added, it becomes the output path
    if(length(benchmarks_out) > 1){
      benchmarks_out <- benchmarks_out[[length(benchmarks_out)]]
      warning("Multiple output tables in output folder, returning newest")
    }
  } else{ # If no new table was added, for out_dir, do nothing. For temp dir, delete. Then throw error.
    if(is.null(out_dir)){
      unlink(outdir, recursive = TRUE)
    }
    stop("No output produced. Check parameters, file paths and input tables.")
  }
  
  message(paste0("Returning BUILDER table: ", benchmarks_out))
  
  out_tab <- utils::read.csv(file.path(outdir, benchmarks_out))
  
  # delete temp directory. If out_dir was provided, don't delete
  if(is.null(out_dir)){
    unlink(outdir, recursive = TRUE)
  }
  
  # make sure BUILDER made some conservation areas
  if(nrow(out_tab) == 0){
    warning("No conservation areas to return, try changing the BUILDER parameters by decreasing intactness and/or area requirements")
    return(NULL)
  }
  
  # return
  out_tab$OID <- NULL
  return(dplyr::as_tibble(out_tab))
}
