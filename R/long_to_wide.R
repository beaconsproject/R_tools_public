### long_to_wide ###
#' Convert long table into wide table with NAs (i.e. BUILDER style table)
#' 
#' Used to prepare long tables for use in BUILDER or functions that take wide tables as 
#' input (e.g.\code{dissolve_catchments_from_table()}).
#' 
#' @param long_df A long style table with values in one column and group names in another.
#' @param col_names The column in long_df from which to get new column names.
#' @param values_col The column in long_df from which to get row values.
#'   
#' @return A wide tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' long_to_wide(
#'   long_df = wide_to_long(builder_table_sample, "group"), 
#'   col_names = "group", 
#'   values_col = "value")
long_to_wide <- function(long_df, col_names, values_col){
  
  # get out table nrow (i.e. longest list of values)
  tbl_rows <- long_df %>%
    dplyr::group_by(.data[[col_names]]) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::summarise(m = max(.data$n)) %>%
    dplyr::pull(.data$m)
  
  values_list <- lapply(unique(long_df[[col_names]]), function(x){
    vals <- long_df[[values_col]][long_df[[col_names]]==x]
    c(vals, rep(NA, tbl_rows - length(vals)))
  })
  names(values_list) <- unique(long_df[[col_names]])
  out_tab <- as.data.frame(do.call(cbind, values_list))
  out_tab <- dplyr::as_tibble(do.call(cbind, values_list))
  
  return(out_tab)
}
