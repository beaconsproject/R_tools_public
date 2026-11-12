### avg_rank ###
#
#' Rank multiple columns and calculate the average rank for each feature.
#'
#' For a given table and a list of columns, calculates the rank of each reserve in each column, and returns the average.
#' 
#' NA values are ranked last. Ties receive the same value, calculated as the average of their collective indices.
#'
#' @param conservation_areas_df table used for ranking.
#' @param columns_asc vector of column names that should be ranked in ascending order
#' @param columns_desc vector of column names that should be ranked in descending order
#'
#' @return A vector of values containing the summed ranks across all requested columns.
#' 
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' tbl <- data.frame(
#'           col1 = c(1,2,3,4,5), 
#'           col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), 
#'           col3 = c(10, NA, 8, 11, NA))
#' avg_rank(tbl, columns_asc = c('col1', 'col3'), columns_desc = 'col2')

avg_rank <- function(conservation_areas_df, columns_asc = c(), columns_desc = c()){
  
  # check columns exist
  if(!all(c(columns_asc, columns_desc) %in% colnames(conservation_areas_df))){
    missing_cols <- c(columns_asc, columns_desc)[!c(columns_asc, columns_desc) %in% colnames(conservation_areas_df)]
    stop(paste0("Columns not in table: ", paste(missing_cols, collapse= ",")))
  }
  
  # NAs - set to max or min value.
  for(column in c(columns_asc, columns_desc)){
    
    if(all(is.na(conservation_areas_df[[column]]))){
      stop(paste0("All values in column ", column, " are NA. Add some values or remove it."))
    }
    if(any(is.na(conservation_areas_df[[column]]))){
      warning(paste0("Column ", column, " contains NA values, they will be ranked last"))
      
      # Set NA to max value plus 1 for asc, min value minus 1 for desc
      if(column %in% columns_asc){
        conservation_areas_df[[column]][is.na(conservation_areas_df[[column]])] <- max(conservation_areas_df[[column]], na.rm = TRUE) + 1
      } 
      if(column %in% columns_desc){
        conservation_areas_df[[column]][is.na(conservation_areas_df[[column]])] <- min(conservation_areas_df[[column]], na.rm = TRUE) - 1
      }
    }
  }
  
  ranks_asc <- lapply(columns_asc, function(x){
    rank(conservation_areas_df[[x]], ties.method = 'average')
  })
  ranks_desc <- lapply(columns_desc, function(x){
    rank(-conservation_areas_df[[x]])
  })
  
  sums <- Reduce("+", c(ranks_asc, ranks_desc))
  count_cols <- length(c(columns_asc, columns_desc))
  return(round(sums / count_cols,  2))
  
}