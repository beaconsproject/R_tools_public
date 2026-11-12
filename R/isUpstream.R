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
