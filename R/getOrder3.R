getOrder3 <- function(catchment_df, current_catchment){
  # get ORDER3 value from the provided CATCHNUM in a data.table
  as.character(catchment_df$ORDER3[catchment_df$CATCHNUM==current_catchment])
}
