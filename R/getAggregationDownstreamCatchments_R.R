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
