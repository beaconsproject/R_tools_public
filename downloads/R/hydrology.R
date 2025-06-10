getOrder1 <- function(catchment_df, current_catchment){
  # get ORDER1 value from the provided CATCHNUM in a data.table
  as.character(catchment_df$ORDER1[catchment_df$CATCHNUM==current_catchment])
}

getOrder2 <- function(catchment_df, current_catchment){
  # get ORDER2 value from the provided CATCHNUM in a data.table
  as.numeric(as.character(catchment_df$ORDER2[catchment_df$CATCHNUM==current_catchment]))
}

getOrder3 <- function(catchment_df, current_catchment){
  # get ORDER3 value from the provided CATCHNUM in a data.table
  as.character(catchment_df$ORDER3[catchment_df$CATCHNUM==current_catchment])
}

getBasin <- function(catchment_df, current_catchment){
  # get BASIN value from the provided CATCHNUM in a data.table
  as.numeric(as.character(catchment_df$BASIN[catchment_df$CATCHNUM==current_catchment]))
}

order3ToOrder2 <- function(order3){
  # replicating builder function to get order 2 from 3
  base64Codes <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "=", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "_", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
  idx0 <- match(substr(order3, 1, 1), base64Codes) - 1 # C# indexes start at zero so minus 1
  idx1 <- match(substr(order3, 2, 2), base64Codes) - 1
  return(idx0 * length(base64Codes) + idx1)
}

utf32_greaterthan <- function(x, y){
  # is x greater than y?
  # used to identify upstream catchments based on comparison of ORDER1 and ORDER1_3 strings.
  
  # get vectors of code points to compare
  xx <- utf8ToInt(x)
  yy <- utf8ToInt(y)
  
  # what is the length of the shortest vector?
  ii <- min(c(length(xx), length(yy)))
  
  # for each index...
  for(i in 1:ii){
    
    # if the value is not the same...
    if(xx[i] != yy[i]){
      
      # test if x > y
      if(xx[i] > yy[i]){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
    # otherwise move on to next index
  }
  
  # If they are the same up to length ii, the upstream catchment will have more characters, so if x is longer than y, return true.
  if(length(xx) > length(yy)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

utf32_greaterthan_vectorized <- function(x, y){
  sapply(x, function(xx){
    utf32_greaterthan(xx, y)
  })
}

isUpstream <- function(catchment_df, current_catchment, other_catchment){
  
  # Test if other_catchment is upstream of current_catchment
  
  # First test: basin = BASIN and order1 = ORDER1 and order2 >= ORDER2
  if(getBasin(catchment_df, other_catchment) == getBasin(catchment_df, current_catchment)){
    if(getOrder1(catchment_df, other_catchment) == getOrder1(catchment_df, current_catchment)){
      if(getOrder2(catchment_df, other_catchment) >= getOrder2(catchment_df, current_catchment)){
        return(TRUE)
      }
    }
  }
  # Second test: basin = BASIN and order1 contains ORDER1 and order1 > ORDER1.ORDER3
  if(getBasin(catchment_df, other_catchment) == getBasin(catchment_df, current_catchment)){
    if(grepl(getOrder1(catchment_df, current_catchment), getOrder1(catchment_df, other_catchment))){
      if(
        utf32_greaterthan(
          getOrder1(catchment_df, other_catchment),
          paste0(getOrder1(catchment_df, current_catchment), getOrder3(catchment_df, current_catchment))
        )
      ){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

getAggregationUpstreamCatchments_BUILDER_method <- function(catchment_df, agg_catchments, neighbours_tib){
  # Get all upstream catchments of an aggregation of catchments (e.g. a list of catchnums in a PA).
  
  # This version copied from BUILDER code.
  # It moves along a stream segment testing only the neighbours of each upstream catchment.
  # If a catchment is upstream it's neighbours are also tested. This has the problem that if the catchments
  # have errors, the upstream progress can be blocked. For example the dissolveing of catchments during construction
  # to eliminate very small catchments can lead to catchments on a stream that should be adjacent, actually not being
  # neighbours.
  # This method also very slow in R.
  # Both problems fixed using the more R-friendly version of this function.
  
  # create inspected list to track upstream catchments that have been inspected
  lst_inspected <- c()
  
  # create list to hold upstream catchments
  lst_upstream <- c()
  
  for(aggCatchmntIdx in agg_catchments){ # for each catchment in the aggregation
    
    # create a pending queue
    pendingQueue <- c()
    
    # add aggCatchmntIdx to the queue
    pendingQueue <- c(pendingQueue, aggCatchmntIdx)
    
    while(length(pendingQueue) > 0){
      
      currentCatchment <- pendingQueue[1]
      
      # remove current catchment from pending queue
      pendingQueue <- pendingQueue[-1]
      
      # get current catchments neighbours
      neighbourCatchments <- neighbours_tib %>%
        dplyr::filter(CATCHNUM == currentCatchment) %>%
        dplyr::pull(neighbours)
      
      for(nbr in neighbourCatchments){
        # if not already inspected
        if(!nbr %in% lst_inspected){
          
          # if not in the aggregation
          if(!nbr %in% agg_catchments){
            
            # test if neighbour is upstream of current agg catchment
            if(isUpstream(catchment_df = catchment_df, current_catchment = aggCatchmntIdx, other_catchment = nbr)){
              
              # if yes, add to upstream list
              lst_upstream <- c(lst_upstream, nbr)
              
              # put this catchment in the queue so its neighbours can be tested
              pendingQueue <- c(pendingQueue, nbr)
              
              # add this catchment to the inspected list
              lst_inspected <- c(lst_inspected, nbr)
            } else{
              # not upstream, do nothing
            }
          } else{
            # if it's in the aggregation, put it in the pending queue so its neighbours can be tested
            #pendingQueue <- c(pendingQueue, nbr)
            #lst_inspected <- c(lst_inspected, nbr)
            
            ### Dropping these lines because if its in the aggregation it'll get tested as an aggCatchmntIdx
            ### The aggCatchmntIdx will still run even if it's already been tested here, so these lines are
            ### unnecessarily repeating the search for neighbors.
          }
        }
      }
    }
  }
  return(lst_upstream)
}

getAggregationUpstreamCatchments_R <- function(catchment_tab, agg_catchments){
  # attempting a more R friendly version with less looping
  
  # This version tests the entire catchment dataset for upstream catchments, repeated
  # for each catchment in the aggregation.
  
  # This is quite slow and many aggregation catchments are repeating the same query because
  # they are on the same stream. A few tricks to speed this up and minimise queries:
  # Remove agg_catchments from the dataframe being queried
  # Remove upstream catchments from the search dataframe so they can't be added multiple times
  # Remove any agg_catchments from the aggregation if they are upstream of a catchment that has already been queried.
  # Start the queries using the most downstream catchments in the aggregation, this will remove as many agg_catchments
  # as possible and minimize the number of queries run.
  
  # make ORDER2 numeric
  catchment_tab$ORDER2 <- as.numeric(as.character(catchment_tab$ORDER2))
  
  # get list of BASINs for the agg_catchments
  basins <- catchment_tab %>%
    dplyr::filter(CATCHNUM %in% agg_catchments) %>%
    dplyr::pull(BASIN) %>%
    unique()
  
  # filter catchments to only include matching BASINs, and to remove agg_catchments
  search_tab <- catchment_tab %>%
    dplyr::filter(!CATCHNUM %in% agg_catchments,
                  BASIN %in% basins)
  
  # make second table to hold just the agg catchments (hopefully faster to query 2 smaller tables than one big)
  agg_tab <- catchment_tab %>%
    dplyr::filter(CATCHNUM %in% agg_catchments) %>%
    dplyr::arrange(nchar(ORDER1), ORDER2) # run the most downstream catchments first to remove as many agg_catchments as possible. This minimizes queries. Rough approximation of downstream is shorter ORDER1 and lower ORDER2 values.
  
  catchList <- c()
  agg_up_list <- c()
  for(current_catchment in agg_tab$CATCHNUM){
    
    #print(paste0(length(catchList), "  -  ", length(agg_up_list), "  -  ", current_catchment))
    
    if(!current_catchment %in% agg_up_list){ # skip if the current agg_catchment is upstream of one that has already been tested. They'll have the same result.
      
      current_BASIN <- getBasin(agg_tab, current_catchment)
      current_ORDER1 <- getOrder1(agg_tab, current_catchment)
      current_ORDER2 <- getOrder2(agg_tab, current_catchment)
      current_ORDER3 <- getOrder3(agg_tab, current_catchment)
      
      # Get all CATCHNUMS matching first test: basin = BASIN and order1 = ORDER1 and order2 >= ORDER2
      test1 <- search_tab %>%
        dplyr::filter(BASIN == current_BASIN & 
                        ORDER1 == current_ORDER1 &
                        ORDER2 >= current_ORDER2) %>%
        dplyr::pull(CATCHNUM)
      
      # Get all CATCHNUMS matching second test: basin = BASIN and order1 contains ORDER1 and order1 > ORDER1.ORDER3
      test2 <- search_tab %>%
        dplyr::filter(BASIN == current_BASIN & 
                        grepl(current_ORDER1, ORDER1) &
                        utf32_greaterthan_vectorized(ORDER1, paste0(current_ORDER1, current_ORDER3))) %>%
        dplyr::pull(CATCHNUM)
      
      # add new upstream catchments to out list
      catchList <- c(catchList, test1, test2)
      
      # remove test1 and test2 catchments from the search table so they can't be added again
      search_tab <- search_tab %>%
        dplyr::filter(!CATCHNUM %in% c(test1, test2))
      
      ###################
      
      # add agg_catchments upstream of current_catchment to a list.
      # these do not need to be tested because they'll have the same result as the current catchment
      test1_agg <- agg_tab %>%
        dplyr::filter(BASIN == current_BASIN &
                        ORDER1 == current_ORDER1 &
                        ORDER2 >= current_ORDER2) %>%
        dplyr::pull(CATCHNUM)
      
      test2_agg <- agg_tab %>%
        dplyr::filter(BASIN == current_BASIN &
                        grepl(current_ORDER1, ORDER1) &
                        utf32_greaterthan_vectorized(ORDER1, paste0(current_ORDER1, current_ORDER3))) %>%
        dplyr::pull(CATCHNUM)
      
      agg_up_list <- c(agg_up_list, test1_agg, test2_agg)
      
    }
  }
  
  # Should already be unique
  outVals <- unique(catchList)
  return(outVals)
}

### get_upstream_catchments ###
#
#' Calculate upstream catchments for a set of input polygons.
#' 
#' Calculates all upstream catchments for each provided protected area polygon and returns as a table, with column names as the unique id 
#' of the protected areas.
#'
#' @param CAs_sf sf object of protected area polygons.
#' @param CAs_id String matching the unique identifier column in \code{CAs_sf}.
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM .
#'
#' @return Tibble where each column name is a unique protected area id, and each row is a catchment making up the 
#' upstream area for that protected area. Blank rows are filled with NA. CATCHNUMs are returned as integers.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' reserves <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
#' get_upstream_catchments(reserves, "network", catchments_sample)
get_upstream_catchments <- function(CAs_sf, CAs_id, catchments_sf){
  
  if(!all(c("ORDER1", "ORDER2", "ORDER3", "BASIN", "CATCHNUM") %in% colnames(catchments_sf))){
    stop("catchments_sf must have attributes: ORDER1, ORDER2, ORDER3, BASIN, CATCHNUM")
  }
  
  # make sure catchnums are integer. Out tables in BUILDER wide format should use integer class for CATCHNUM
  catchments_sf <- make_catchnum_integer(catchments_sf)
  
  # get list of catchnums in each PA
  CAs_catchnums_tab <- catchnums_in_polygon(CAs_sf, CAs_id, catchments_sf)
  
  up_agg_list <- list()
  for(col_id in colnames(CAs_catchnums_tab)){
    
    # get list of catchments
    agg_catchments <- get_catch_list(col_id, CAs_catchnums_tab)
    up_agg <- getAggregationUpstreamCatchments_R(catchments_sf, agg_catchments)
    
    # add to up_agg_list
    up_agg_list[[col_id]] <- up_agg
  }
  
  out_df <- dplyr::as_tibble(list_to_wide(up_agg_list))
  
  return(out_df)
}


getAggregationDownstreamCatchments_R <- function(catchment_tab, agg_catchments){
  # attempting a more R friendly version with less looping
  
  # ORDER 1 of a new upstream segment at a stream fork is concatenation of [ORDER1, ORDER3, '1'] of the last downstream segment
  # i.e. for and ORDER1 stream of .003U1, the next downstream segment on the main river stem will have ORDER 1 and 3 of '.00' 
  # and '3U'.
  
  # For any ORDER 1 segment, we can therefore find all combinations of ORDER1 ORDER3 segments at the downstream intersections,
  # i.e the 'next downstream' catchments.
  # In the example of '.003U10e1' there would be two 'next downstream' catchments with ORDER values: '.003U1' + '0e' and '.00' + '3U'.
  
  # For the original segment, and each of these 'next downstream' segments, we can use the ORDER2 value from the catchments table to
  # grab all downstream segments on that ORDER1 stream.
  
  # Note that because many catchments are dissolved, the exact 'next downstream' ORDER3 catchment on a given ORDER1 stream might not exist.
  # Instead of searching for that specific catchment, we convert the ORDER3 to it's matching ORDER2 value and simply query on all lower
  # ORDER2 values. It doesn't matter that the actual ORDER2 catchment doesn't exist, we just want all smaller values.
  
  # This is quite slow and many aggregation catchments are repeating the same query because
  # they are on the same stream. A few tricks to speed this up and minimise queries:
  # Remove agg_catchments from the dataframe being queried
  # Remove downstream catchments from the search dataframe so they can't be added multiple times
  # Remove any agg_catchments from the aggregation if they are downstream of a catchment that has already been queried.
  # Start the queries using the most upstream catchments in the aggregation, this will remove as many agg_catchments
  # as possible and minimize the number of queries run.
  
  # make ORDER2 numeric
  catchment_tab$ORDER2 <- as.numeric(as.character(catchment_tab$ORDER2))
  
  # filter catchments to remove agg_catchments
  search_tab <- catchment_tab %>%
    dplyr::filter(!CATCHNUM %in% agg_catchments)
  
  # make second table to hold just the agg catchments (hopefully faster to query 2 smaller tables than one big)
  agg_tab <- catchment_tab %>%
    dplyr::filter(CATCHNUM %in% agg_catchments) %>%
    dplyr::arrange(dplyr::desc(nchar(ORDER1)), dplyr::desc(ORDER2)) # run the most upstream catchments first to remove as many agg_catchments as possible. This minimizes queries. Rough approximation of downstream is longest ORDER1 and highest ORDER2 values.
  
  catchList <- c()
  for(current_catchment in agg_tab$CATCHNUM){
    
    if(!current_catchment %in% catchList){ # skip if the current agg_catchment is downstream of one that has already been tested. They'll have the same result.
      
      # Get all downstream on current catchments ORDER1
      current_ORDER1 <- getOrder1(agg_tab, current_catchment)
      current_ORDER2 <- getOrder2(agg_tab, current_catchment)
      current_BASIN <- getBasin(agg_tab, current_catchment)
      
      current_catchList <- search_tab %>%
        dplyr::filter(BASIN == current_BASIN &
                        ORDER1 == current_ORDER1 &
                        ORDER2 < current_ORDER2) %>%
        dplyr::pull(CATCHNUM)
      
      # Get ORDER1 and ORDER2 for each 'next downstream' catchments
      # Use them to grab all lower ORDER2 catchments on the stream and add to catchList
      count <- (nchar(current_ORDER1)/3) - 1
      for(i in 1:count){
        next_ORDER1 <- substr(current_ORDER1, 1, nchar(current_ORDER1) - i*3)
        next_ORDER3 <- substr(current_ORDER1, (nchar(current_ORDER1) - i*3)+1, (nchar(current_ORDER1) - i*3) + 2)
        next_ORDER2 <- order3ToOrder2(next_ORDER3)
        
        down_i <- search_tab %>%
          dplyr::filter(BASIN == current_BASIN &
                          ORDER1 == next_ORDER1 &
                          ORDER2 <= next_ORDER2) %>%
          dplyr::pull(CATCHNUM)
        
        current_catchList <- c(current_catchList, down_i)
      }
      
      # add new downstream catchments to out list
      catchList <- c(catchList, current_catchList)
      
      # remove new downstream catchments from the search table so they can't be added again
      if(length(current_catchList) > 0){
        search_tab <- search_tab %>%
          dplyr::filter(!CATCHNUM %in% current_catchList)
      }
    }
  }
  
  # Should already be unique
  outVals <- unique(catchList)
  return(outVals)
}

### get_downstream_catchments ###
#
#' Calculate downstream catchments for a set of input polygons.
#' 
#' Calculates all downstream catchments for each provided protected area polygon and returns as a table, with column names as the unique id 
#' of the protected areas.
#'
#' @param CAs_sf sf object of protected area polygons.
#' @param CAs_id String matching the unique identifier column in \code{CAs_sf}.
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM .
#'
#' @return Tibble where each column name is a unique protected area id, and each row is a catchment making up the 
#' upstream area for that protected area. Blank rows are filled with NA. CATCHNUMs are returned as integers.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' reserves <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
#' get_downstream_catchments(reserves, "network", catchments_sample)
get_downstream_catchments <- function(CAs_sf, CAs_id, catchments_sf){
  
  if(!all(c("ORDER1", "ORDER2", "ORDER3", "BASIN", "CATCHNUM") %in% colnames(catchments_sf))){
    stop("catchments_sf must have attributes: ORDER1, ORDER2, ORDER3, BASIN, CATCHNUM")
  }
  
  # make sure catchnums are integer. Out tables in BUILDER wide format should use integer class for CATCHNUM
  catchments_sf <- make_catchnum_integer(catchments_sf)
  
  # get list of catchnums in each PA
  CAs_catchnums_tab <- catchnums_in_polygon(CAs_sf, CAs_id, catchments_sf)
  
  down_agg_list <- list()
  for(col_id in colnames(CAs_catchnums_tab)){
    
    # get list of catchments
    agg_catchments <- get_catch_list(col_id, CAs_catchnums_tab)
    down_agg <- getAggregationDownstreamCatchments_R(catchments_sf, agg_catchments)
    
    # add to up_agg_list
    down_agg_list[[col_id]] <- down_agg
  }
  
  out_df <- dplyr::as_tibble(list_to_wide(down_agg_list))
  
  return(out_df)
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
#' @param CAs_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param CAs_id character specify unique identifier. Default is network.
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
calc_dci <- function(CAs_sf, stream_sf, CAs_id = "network", buffer_width = 0.1){
  
  stopifnot(sf::st_crs(CAs_sf) == sf::st_crs(stream_sf))
  
  # set stream and conservation area attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(CAs_sf) = "constant"
  
  # clip streams to full area of all conservation areas
  reserve_dci <- CAs_sf %>%
    sf::st_intersection(stream_sf, sf::st_union()) %>% # get streams just for the required region
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::summarise(geometry = sf::st_union(geometry)) %>% # merge into single feature
    sf::st_intersection(CAs_sf) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select(dplyr::all_of(CAs_id)) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(CAs_id)) %>% # for each network...
    dplyr::summarise(L = sum(stream_length), dci = sum((stream_length*stream_length) / (L*L)))  # calculate L2 then use to calculate dci
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  dci <- dplyr::left_join(sf::st_drop_geometry(CAs_sf), reserve_dci, by = CAs_id) %>%
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
#' @param CAs_sf sf object of conservation areas in which to calculate DCI.
#' @param stream_sf sf object of river network. Must have streams grouped in a BASIN attribute.
#' @param CAs_id character specify unique identifier. Default is network.
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

calc_lwdci <- function(CAs_sf, stream_sf, CAs_id = "network", buffer_width = 0.1){
  
  stopifnot(sf::st_crs(CAs_sf) == sf::st_crs(stream_sf))
  # set stream and reserve attributes to be constant throughout each polygon to avoid warnings in st_intersection
  sf::st_agr(stream_sf) = "constant"
  sf::st_agr(CAs_sf) = "constant"
  
  # clip streams to full area of all reserves
  stream_prepped <- 
    sf::st_intersection(stream_sf, sf::st_union(CAs_sf$geometry)) %>%
    dplyr::filter(BASIN != -1) %>% # remove isolated stream segments
    sf::st_buffer(dist = buffer_width, endCapStyle = "ROUND") %>% # buffer to make sure streams are connected
    dplyr::group_by(BASIN) %>%
    dplyr::summarise(geometry = sf::st_union(geometry)) # merge into single feature per BASIN
  
  sf::st_agr(stream_prepped) = "constant"
  
  reserve_lwdci <- sf::st_intersection(CAs_sf, stream_prepped) %>% # intersect with reserves to get buffered stream for each reserve
    dplyr::select(dplyr::all_of(CAs_id), BASIN) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% # this is needed to avoid geometries being lost in the POLYGON cast
    sf::st_cast("POLYGON", warn = FALSE) %>% # explode into individual stream segments
    dplyr::mutate(stream_length = as.numeric(sf::st_area(geometry)) / buffer_width) %>% # divide area by buffer to get length of each stream segment
    sf::st_drop_geometry() %>% # drop the geometry for speed
    dplyr::group_by(!!sym(CAs_id), BASIN) %>% # for each network and basin...
    dplyr::summarise(L_basin = sum(stream_length), 
                     dci_basin = sum((stream_length*stream_length) / (L_basin*L_basin)), .groups = "drop") %>% # calculate L2 then use to calculate dci
    dplyr::group_by(!!sym(CAs_id)) %>%
    dplyr::summarise(lwdci = sum((L_basin / sum(L_basin)) * dci_basin), .groups = "drop") # calc lwdci. Weight is proportion of reserve streams in each basin
  
  # reserves that do not intersect the stream network get dropped during st_intersection.
  # join dci back to original reserves and set missing reserves to have dci of 0
  lwdci <- dplyr::left_join(sf::st_drop_geometry(CAs_sf), reserve_lwdci, by = CAs_id) %>%
    tidyr::replace_na(list(lwdci=0)) %>%
    dplyr::pull(lwdci) %>%
    round(3)
  
  return(lwdci)
}
