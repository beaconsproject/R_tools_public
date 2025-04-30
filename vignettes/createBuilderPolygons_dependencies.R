dissolve_catchments_from_table <- function(catchments_sf, input_table, out_feature_id=NULL, calc_area = FALSE, intactness_id = NULL, dissolve_list = c(), drop_table = NULL){
  check_catchnum(catchments_sf) # check for CATCHNUM
  check_for_geometry(catchments_sf)
  #check_catchnum_class(catchments_sf, input_table) # Check catchments match, warning if not
  input_table <- remove_oid(input_table) #drop oid column is it exists
  check_for_rows(input_table)
  
  # get colnames to process
  if(length(dissolve_list > 0)){
    feature_list <- dissolve_list
  } else{
    feature_list <- colnames(input_table)
  }
  
  saveCount <- 1
  for(col_id in feature_list){
    
    # separate names if combination
    col_ids <- sep_network_names(col_id)
    
    # get list of catchments
    catchments_list <- get_catch_list(col_ids, input_table)
    
    # only proceed if there are catchments in the list. if no catchments, the network will not be included in the output table
    if(length(catchments_list) > 0){
      
      # drop catchments if requested
      if(!is.null(drop_table)){
        drop_list <- get_catch_list(col_ids, drop_table)
        
        catchments_list <- catchments_list[!catchments_list %in% drop_list]
      }
      
      # dissolve based on parameters
      if(calc_area){
        dslv <- catchments_sf %>%
          dplyr::filter(CATCHNUM %in% catchments_list) %>%
          dplyr::summarise(geometry = sf::st_union(geometry)) %>%
          sf::st_buffer(20) %>%
          sf::st_buffer(-20) %>%
          dplyr::mutate(id = col_id,
                        area_km2 = round(as.numeric(sf::st_area(geometry) / 1000000), 2))
      } else{
        dslv <- catchments_sf %>%
          dplyr::filter(CATCHNUM %in% catchments_list) %>%
          dplyr::summarise(geometry = sf::st_union(geometry)) %>%
          sf::st_buffer(20) %>%
          sf::st_buffer(-20) %>%
          dplyr::mutate(id = col_id)
      }
      
      # join AWI if requested
      if(!is.null(intactness_id)){
        if(intactness_id %in% colnames(catchments_sf)){
          awi <- catchments_sf %>%
            dplyr::filter(CATCHNUM %in% catchments_list) %>%
            dplyr::mutate(area = as.numeric(sf::st_area(geometry))) %>%
            sf::st_drop_geometry() %>%
            dplyr::summarise(AWI = sum(.[[intactness_id]] * area) / sum(area))
          
          dslv$AWI <- round(awi$AWI, 4)
        } else{
          warning(paste0("Area-weighted intactness cannot be calculated, ", '"', intactness_id, '"', " not in catchments_sf"))
        }
      }
      
      # set out name
      if(!is.null(out_feature_id)){
        names(dslv)[names(dslv) == "id"] <- out_feature_id
      }
      
      # reorder columns - move geometry to last
      dslv <- dslv  %>%
        dplyr::relocate(geometry, .after = dplyr::last_col())
      
      # append to df
      if(saveCount == 1){
        out_sf <- dslv
        saveCount <- saveCount + 1
      } else{
        out_sf <- rbind(out_sf, dslv)
      }
    }
  }
  return(out_sf)
}
calc_dci <- function(conservation_area_sf, stream_sf, col_name = "network", buffer_width = 0.1){
  
  stopifnot(sf::st_crs(conservation_area_sf) == sf::st_crs(stream_sf))
  
  # set stream and conservation area attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(conservation_area_sf) = "constant"
  
  # clip streams to full area of all conservation areas
  reserve_dci <- conservation_area_sf %>%
    sf::st_intersection(stream_sf, sf::st_union()) %>% # get streams just for the required region
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::summarise(geometry = sf::st_union(geometry)) %>% # merge into single feature
    sf::st_intersection(conservation_area_sf) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select(dplyr::all_of(col_name)) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(col_name)) %>% # for each network...
    dplyr::summarise(L = sum(stream_length), dci = sum((stream_length*stream_length) / (L*L)))  # calculate L2 then use to calculate dci
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  dci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_dci, by = col_name) %>%
    tidyr::replace_na(list(dci=0)) %>%
    dplyr::pull(dci) %>%
    round(3)
  
  return(dci)
}
calc_lwdci <- function(conservation_area_sf, stream_sf, col_name = "network", buffer_width = 0.1){
  
  stopifnot(sf::st_crs(conservation_area_sf) == sf::st_crs(stream_sf))
  # set stream and reserve attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(conservation_area_sf) = "constant"
  
  # clip streams to full area of all reserves
  stream_prepped <- 
    sf::st_intersection(stream_sf, sf::st_union(conservation_area_sf$geometry)) %>%
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::group_by(BASIN) %>%
    dplyr::summarise(geometry = sf::st_union(geometry)) # merge into single feature per BASIN
  
  sf::st_agr(stream_prepped) = "constant"
  
  reserve_lwdci <- sf::st_intersection(conservation_area_sf, stream_prepped) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select(dplyr::all_of(col_name), BASIN) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(col_name), BASIN) %>% # for each network and basin...
    dplyr::summarise(L_basin = sum(stream_length), 
                     dci_basin = sum((stream_length*stream_length) / (L_basin*L_basin)), .groups = "drop") %>% # calculate L2 then use to calculate dci
    dplyr::group_by(!!sym(col_name)) %>%
    dplyr::summarise(lwdci = sum((L_basin / sum(L_basin)) * dci_basin), .groups = "drop") # calc lwdci. Weight is proportion of reserve streams in each basin
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  lwdci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_lwdci, by = col_name) %>%
    tidyr::replace_na(list(lwdci=0)) %>%
    dplyr::pull(lwdci) %>%
    round(3)
  
  return(lwdci)
}
sep_network_names <- function(network_names){
  
  if(length(network_names) == 1){
    out_val <- strsplit(network_names, "__")[[1]]
  }
  
  if(length(network_names) > 1){
    out_val <- lapply(network_names, function(x){
      strsplit(x, "__")[[1]]
    })
    names(out_val) <- network_names
  }
  return(out_val)
}
check_catchnum <- function(catchments_sf){
  
  # check CATCHNUM exists
  if(!"CATCHNUM" %in% names(catchments_sf)){
    stop("Catchments must contain column 'CATCHNUM'")
  }
}
check_catchnum_class <- function(catchments_sf, builder_table){
  col_classes <- sapply(colnames(builder_table), function(x) class(builder_table[[x]]))
  if(!all(col_classes == class(catchments_sf$CATCHNUM))){
    warning(paste0("Table column classes do not match class(catchments_sf$CATCHNUM) for columns: ", paste0(colnames(builder_table)[col_classes != class(catchments_sf$CATCHNUM)], collapse=", ")))
  }
}
get_catch_list <- function(pa_ids, input_table){
  
  # check nets are in input_table
  if(!all(pa_ids %in% colnames(input_table))){
    stop(paste0("names do not appear as colnames in the input table: ", pa_ids[!pa_ids %in% colnames(input_table)]))
  }
  
  # get unique catchments vector from input_table using vector of colnames
  net_catchments <- unique(unlist(lapply(pa_ids, function(x){
    input_table[[x]]
  })))
  
  net_catchments[!is.na(net_catchments)]
}
check_for_geometry <- function(in_sf){
  
  if(!"geometry" %in% names(in_sf)){
    stop("Object must contain column: geometry")
  }
}
remove_oid <- function(in_table){
  if("OID" %in% colnames(in_table)){
    out_table <- in_table %>%
      dplyr::select(-OID)
  } else{
    out_table <- in_table
  }
  return(out_table)
}
check_for_rows <- function(in_table){
  if(nrow(in_table) == 0){
    stop("input_table has no data")
  }
}

prep_input_column <- function(out_dir, type= "BENCHMARKS"){
  if(is.null(out_dir)){
    stop("You need to point to a Builder_output directory")
  }
  
  if(type =="BENCHMARKS"){
    benchmarks_out <- list.files(out_dir, pattern = "COLUMN_All_Unique_BAs.csv")
  }else if (type == "UPSTREAM"){
    benchmarks_out <- list.files(out_dir, pattern = "UPSTREAM_CATCHMENTS_COLUMN.csv")
  }else if (type == "DOWNSTREAM"){
    benchmarks_out <- list.files(out_dir, pattern = "DOWNSTREAM_CATCHMENTS_COLUMN.csv")
  }else{
    stop("Wrong Type. Please select BENCHMARK, UPSTREAM or DOWNSTREAM")
  }
  
  if(length(benchmarks_out) > 1){
    benchmarks_out <- benchmarks_out[[length(benchmarks_out)]]
    warning("Multiple output tables in output folder, returning newest")
  } else if (length(benchmarks_out)==0){ # If no new table was added, for out_dir, do nothing. For temp dir, delete. Then throw error.
    stop("No output produced. Check parameters, file paths and input tables.")
  }
  
  message(paste0("Returning BUILDER table: ", benchmarks_out))
  out_tab <- utils::read.csv(file.path(out_dir, benchmarks_out))
  
  # make sure BUILDER made some conservation areas
  if(nrow(out_tab) == 0){
    warning("No conservation areas to return, try changing the BUILDER parameters by decreasing intactness and/or area requirements")
    return(NULL)
  }
  
  # return
  out_tab$OID <- NULL
  return(dplyr::as_tibble(out_tab))
}
