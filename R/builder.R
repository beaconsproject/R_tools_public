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
#' @param conservation_areas_intactness Minimum area-weighted intactness of final conservation areas. Only conservation areas meeting this value will be returned.
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
                    catchment_level_intactness = 1, conservation_areas_intactness = 1, area_target_proportion = 1, # main parameters
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
    missing_cols <- check_colnames(seeds, c("CATCHNUM", "Areatarget"))
  }else if(data_source == "Reserve"){
    missing_cols <- check_colnames(seeds, c(reserve_name, "Areatarget", "CATCHNUM"))
  }else{
    stop("data_source must be either catchment or Reserve")
  }
  if (length(missing_cols) > 0) stop("Column(s) missing from seeds: ", paste(missing_cols, collapse = ", "))

  check_seeds_areatargets(seeds)
  check_seeds_in_catchments(seeds, catchments_sf)
  utils::write.csv(seeds, file.path(outdir, "seeds.txt"), row.names = FALSE)
  
  # save neighbours table to outdir
  message("checking neighbours table")
  missing_cols <- check_colnames(neighbours, c("CATCHNUM", "neighbours", "key"))
  if (length(missing_cols) > 0) stop("Column(s) missing from neighbours: ", paste(missing_cols, collapse = ", "))
  neighbours <- make_all_integer(neighbours)
  utils::write.csv(neighbours, file.path(outdir, "neighbours.csv"), row.names = FALSE)
  
  # save catchment attributes
  catchments_csv <- sf::st_drop_geometry(catchments_sf)
  missing_cols <- check_colnames(catchments_csv, c(area_land, area_water, skeluid, catchnum, subzone, zone, basin, order1, order2, order3, stream_length, intactness, isolated)) # check all provided columns are present
  if (length(missing_cols) > 0) stop("Column(s) missing from catchments_sf: ", paste(missing_cols, collapse = ", "))
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
                          conservation_areas_intactness,
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
