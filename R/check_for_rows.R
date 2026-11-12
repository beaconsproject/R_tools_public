# Check for rows
check_for_rows <- function(in_table){
  if(nrow(in_table) == 0){
    stop("input_table has no data")
  }
}
