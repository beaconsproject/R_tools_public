getOrder2 <- function(catchment_df, current_catchment){
  # get ORDER2 value from the provided CATCHNUM in a data.table
  as.numeric(as.character(catchment_df$ORDER2[catchment_df$CATCHNUM==current_catchment]))
}
