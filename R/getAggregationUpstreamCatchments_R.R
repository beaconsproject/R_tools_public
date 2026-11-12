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
