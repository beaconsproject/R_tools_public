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
#' @importFrom rlang .data
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

genReserveSeed <- function (catchments, reserves, reserveName, areaTarget, joinType, outDir){
  options(scipen = 999)
  #browser()
  # JOIN TYPE
  if(joinType =="INTERSECT"){
    catch_join <- st_intersects(reserves, catchments)
  }else if(joinType =="CENTROID"){
    #browser()
    catch_centroids <- st_centroid(catchments)
    catch_join <- st_intersects(reserves, catch_centroids, sparse = TRUE)
  }
  
  out_tab <- data.frame(
    reservename = reserves[[reserveName]],
    areatarget = areaTarget,
    catchnums  = sapply(catch_join, function(idx) {
      paste(catchments$CATCHNUM[idx], collapse = ",")
    })
  ) 
  
  # # Split the 'catchnums' column into a list of vectors
  catchnum_list <- strsplit(out_tab$catchnums, ",")
  
  # Determine the maximum number of catch numbers
  max_catchnums <- max(sapply(catchnum_list, length))
  
  # Create a row with the specified values
  row <- c("reservename", "areatarget", "catchnums", rep(NA, max_catchnums-1))
  
  # Create the data frame
  seed_init <- as.data.frame(t(row), stringsAsFactors = FALSE)
  utils::write.table(seed_init, file = file.path(outDir, "Builder_input/seeds.csv"), sep= ",", na= "", row.names = FALSE, col.names= FALSE)
  
  # Specify the output file path
  output_file <- file.path(outDir, "Builder_input/seeds.csv")
  
  # Open a connection to the output file
  file_conn <- file(output_file, open = "a")
  
  # Pad shorter vectors with NA
  catchnum_list <- lapply(catchnum_list, function(x) {
    length(x) <- max_catchnums
    return(x)
  })
  
  # Iterate over each row in the data list
  for (i in 1:length(catchnum_list)) {
    reserveName <- as.character(out_tab[i, 1])
    areatarget <- as.character(out_tab[i, 2])
    catchnums <- catchnum_list[[i]]
    # Combine the elements into a single comma-separated string
    last_non_na <- max(which(catchnums != ""))
    catchnums <- catchnums[1:last_non_na]
    
    line <- paste(c(reserveName, areatarget, catchnums), collapse = ",")
    
    # Write the line to the file
    writeLines(line, con = file_conn)
  }
  
  # Close the file connection
  close(file_conn)
  #browser()
  out_tab <- utils::read.csv(output_file)
  return(dplyr::as_tibble(out_tab))
  return(seed)
}

#Generate polygone
gen_poly_from_rowTable <- function(catchments, outDir, type = "UPSTREAM"){
  flist <- list.files(path = file.path(outDir,"Builder_output"), pattern = type)
  last_file <- flist[length(flist)]
  
  bfile <- readr::read_lines(file.path(outDir,"Builder_output",last_file), skip=1)

  # Step 2: Split each line by comma
  bfile_split <- str_split(bfile, ",", simplify = FALSE)
  
  # Step 3: Drop second element of each row
  bfile_clean <- lapply(bfile_split, function(x) x[-2])
  
  # Step 4: Convert each vector to a data frame (to work with rbind.fill)
  bfile_df_list <- lapply(bfile_clean, function(x) as.data.frame(t(x), stringsAsFactors = FALSE))
  
  # Step 5: Use rbind.fill to combine the data frames
  bfile_combined <- plyr::rbind.fill(bfile_df_list)
  
  # Step 6: Transpose and convert to tibble
  bfile_long <- as_tibble(t(bfile_combined), .name_repair = "minimal")
  
 # Set column names from first row
  colnames(bfile_long) <- bfile_long[1, ]
  bfile_long <- bfile_long[-1, ]
  bfile_long <-  bfile_long %>%
    mutate_all(as.numeric)
  
  sf_poly <- suppressWarnings(dissolve_catchments_from_table(catchments_sf, bfile_long))
  sf_poly <- sf_poly %>%
    st_buffer(dist = 20) %>% 
    st_buffer(dist = -20)
  return(sf_poly)
}


dissolve_catchments_from_table <- function(catchments_sf, input_table, out_feature_id=NULL, calc_area = FALSE, intactness_id = NULL, dissolve_list = c(), drop_table = NULL){
  check_catchnum(catchments_sf) # check for CATCHNUM
  check_for_geometry(catchments_sf)
  check_catchnum_class(catchments_sf, input_table) # Check catchments match, warning if not
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
          dplyr::mutate(id = col_id,
                        area_km2 = round(as.numeric(sf::st_area(geometry) / 1000000), 2))
      } else{
        dslv <- catchments_sf %>%
          dplyr::filter(CATCHNUM %in% catchments_list) %>%
          dplyr::summarise(geometry = sf::st_union(geometry)) %>%
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

make_catchnum_integer <- function(catchments_sf){
  
  # check CATCHNUM exists
  if(!"CATCHNUM" %in% names(catchments_sf)){
    stop("Catchments must contain column 'CATCHNUM'")
  }
  
  catchments_sf$CATCHNUM <- as.integer(catchments_sf$CATCHNUM)
  
  return(catchments_sf)
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

check_for_geometry <- function(in_sf){
  
  if(!"geometry" %in% names(in_sf)){
    stop("Must contain column: geometry")
  }
}

# Check for rows
check_for_rows <- function(in_table){
  if(nrow(in_table) == 0){
    stop("input_table has no data")
  }
}

# remove OID from tables comgin out of BUILDER
remove_oid <- function(in_table){
  if("OID" %in% colnames(in_table)){
    out_table <- in_table %>%
      dplyr::select(-OID)
  } else{
    out_table <- in_table
  }
  return(out_table)
}

# sep_network_names
#' Split network names into their constituent parts.
#'
#' Takes a vector of network names in the form output by [gen_network_names()] (i.e. \code{"__"} separated),
#' and splits the network names back into their original names.
#' 
#' If a single name is provided, a vector of names is returned. If a vector with length >1 is provided, a list
#' of vectors is returned with each list element named with the network name.
#'
#' @param network_names String of a network name, or a vector of multiple network names. Must include the \code{"__"} separator.
#'
#' @return A vector of names, or a list of vectors.
#' @export
#'
#' @examples
#' sep_network_names("PB_1__PB_2")
#' sep_network_names(c("PB_1__PB_2", "PB_1__PB_3"))
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2)
#' sep_network_names(x)
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 4)
#' sep_network_names(x)

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


### get_catch_list ###
#
#' Get unique sets of catchment values from a table.
#' 
#' Designed for use with wide format tables such as those output by \code{get_upstream_catchments} where each 
#' column name is a unique protected area id, and each row is a catchment id. Catchment values from all columns 
#' included in pa_ids are combined and returned after removing duplicates.
#' 
#' @param pa_ids A single column name in input_table, or a vector of multiple column names.
#' @param input_table Data frame where column names are protected area names and rows are catchments making 
#' up each polygon. e.g. output from [get_upstream_catchments()].
#'   
#' @return A vector of unique catchment values.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' get_catch_list("PB_0001", builder_table_sample)
#' get_catch_list(c("PB_0001", "PB_0003"), builder_table_sample)
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

### calc_dci ###
#
#' Calculate dendritic connectivity index (DCI) for a set of input polygons.
#' 
#' Values range between 0 and 1, with 1 indicating a completely connected river network.
#'
#' Dendritic connectivity measures the longitudinal connectivity of a river network inside a set of reserves
#' based on the expected probability of an organism being able to move between two random points in the
#' network [(Cote et al 2009)](https://link.springer.com/article/10.1007/s10980-008-9283-y).
#'
#' DCI = sum ( li2 / L2 )
#' li - length of stream section
#' L - total length of all stream sections in conservation area
#' 
#'
#' @param conservation_area_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param buffer_width Width of buffer to apply to stream segments. Defaults to 0.1. Used to
#' ensure adjacent stream segments are connected during analysis.
#'
#' @return Vector of numeric DCI values matching the input features. 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' calc_dci(conservation_areas, streams_sample)
#' 
calc_dci <- function(conservation_area_sf, stream_sf, column_id = "network", buffer_width = 0.1){
  #browser()
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
    dplyr::select({{column_id}}) %>% # drop attributes except network id
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(column_id)) %>% # for each network...
    dplyr::summarise(L = sum(stream_length), dci = sum((stream_length*stream_length) / (L*L)))  # calculate L2 then use to calculate dci
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  dci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_dci, by = column_id) %>%
    tidyr::replace_na(list(dci=0)) %>%
    dplyr::pull(dci) %>%
    round(3)
  
  return(dci)
}

#' calc_lwdci
#' 
#' Calculate stream network length-weighted dendritic connectivity index (DCI) for a set of input polygons.
#'
#' Dendritic connectivity measures the longitudinal connectivity of a river network inside a set of conservation areas
#' based on the expected probability of an organism being able to move between two random points in the
#' network [(Cote et al 2009)](https://link.springer.com/article/10.1007/s10980-008-9283-y).
#'
#' DCI = sum ( li2 / L2 )
#' li - length of stream section
#' L - total length of all stream sections in conservation area
#' 
#' Values range between 0 and 1, with 1 indicating a completely connected river network.
#' 
#' For length weighted DCI, a separate DCI measure is calculated for each group of BASIN streams 
#' within the conservation area, then a weighted average is calculated where the weights are the lengths 
#' of streams in each BASIN.
#'
#' @param conservation_area_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param buffer_width Width of buffer to apply to stream segments. Defaults to 0.1. Used to
#' ensure adjacent stream segments are connected during analysis.
#'
#' @return Vector of numeric DCI values matching the input features. 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' calc_lwdci(conservation_areas, streams_sample)
#' 

calc_lwdci <- function(conservation_area_sf, stream_sf, column_id = "network", buffer_width = 0.1){
  
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
    dplyr::select({{column_id}}, BASIN) %>% # drop attributes except network and BASIN
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(column_id), BASIN) %>% # for each network and basin...
    dplyr::summarise(L_basin = sum(stream_length), 
                     dci_basin = sum((stream_length*stream_length) / (L_basin*L_basin)), .groups = "drop") %>% # calculate L2 then use to calculate dci
    dplyr::group_by(!!sym(column_id)) %>%
    dplyr::summarise(lwdci = sum((L_basin / sum(L_basin)) * dci_basin), .groups = "drop") # calc lwdci. Weight is proportion of reserve streams in each basin
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  lwdci <- dplyr::left_join(sf::st_drop_geometry(conservation_area_sf), reserve_lwdci, by = column_id) %>%
    tidyr::replace_na(list(lwdci=0)) %>%
    dplyr::pull(lwdci) %>%
    round(3)
  
  return(lwdci)
}

