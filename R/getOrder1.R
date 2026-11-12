getOrder1 <- function(catchment_df, current_catchment){
  # get ORDER1 value from the provided CATCHNUM in a data.table
  as.character(catchment_df$ORDER1[catchment_df$CATCHNUM==current_catchment])
}
