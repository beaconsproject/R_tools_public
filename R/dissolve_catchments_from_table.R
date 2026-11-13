### dissolve_catchments_from_table ###
#
#' Get lists of catchments from a table and dissolve.
#' 
#' Takes lists of catchments from a table where each column represents a polygon to be created, and dissolves the 
#' catchments into a single polygon feature. Optionally adds area in km2, and area-weighted intactness.
#' 
#' \code{input_table} is in the form exported from the beaconsbuilder package, where each column lists catchment ids, 
#' and each column name represents the unique identifier of the output polygon.
#' 
#' This function is used to create conservation area polygons using tables created by the beaconsbuilder package, 
#' and upstream polygons using tables created by [get_upstream_catchments()].
#' 
#' The optional \code{dissolve_list} parameter can be used to filter the column names in \code{input_table} that will
#' be dissolved. The \code{dissolve_list} can also contain combinations of column name separated by \code{__}, in this
#' case the output polygon will dissolve the combined area of catchments from all columns in the string. For example, 
#' the string "PA1__PA2" in the \code{dissolve_list} will output a polygon representing the dissolved area of all
#' catchments listed in the PA1 and PA2 columns of \code{input_table}.
#' 
#' The main application of \code{drop_table} is to remove conservation area catchments from the dissolved area when dissolving
#' upstream catchments tables. See vignette('overview') for an example.
#' 
#' @param catchments_sf sf object of the catchments dataset with unique identifier column: CATCHNUM.
#' @param input_table Data frame where column names are polygon names and rows are catchments making up each polygon.
#' e.g. output from [get_upstream_catchments()].
#' @param out_feature_id String representing the output column name holding the polygon unique identifiers.
#' @param calc_area If TRUE, an area_km2 column is added to the output containing the polygon area.
#' @param intactness_id Optional string identifying an intactness column (values between 0 and 1) in catchments_sf. If provided, 
#' used to calculate the area weighted intactness (AWI) of the dissolved polygons.
#' @param dissolve_list Vector of columns in \code{input_table} to include in the output. Defaults to colnames(input_table).
#'   Can also dissolve multiple columns from \code{input_table} together by combining column names with \code{__} (e.g. PB_0001__PB_0002).
#' @param drop_table A table in the same format \code{input_table} listing catchments to be dropped from the dissolve area.
#'   
#' @return A sf object of polygons with unique identifier column \code{out_feature_id}.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library(sf)
#' 
#' # Individual conservation areas
#' conservation_areas <- dissolve_catchments_from_table(
#'   catchments_sf = catchments_sample, 
#'   input_table = builder_table_sample, 
#'   out_feature_id = "network",
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' plot(conservation_areas)
#' 
#' # conservation areas combined into a network
#' network <- dissolve_catchments_from_table(
#'   catchments_sf = catchments_sample, 
#'   input_table = builder_table_sample, 
#'   out_feature_id = "network",
#'   dissolve_list = c("PB_0001__PB_0002"))
#' plot(network)
#' 
#' # Upstream polygons
#' upstream_table <- get_upstream_catchments(conservation_areas, "network", catchments_sample)
#' dissolve_catchments_from_table(
#'   catchments_sample, 
#'   upstream_table, 
#'   "network",
#'   TRUE)
#'   
dissolve_catchments_from_table <- function(catchments_sf, input_table, out_feature_id=NULL, calc_area = FALSE, intactness_id = NULL, dissolve_list = c(), drop_table = NULL){
  missing_cols <- check_colnames(catchments_sf, "CATCHNUM")
  if (length(missing_cols) > 0) stop("Column(s) missing from catchments_sf: ", paste(missing_cols, collapse = ", "))
  check_for_geometry(catchments_sf)
  input_table <- remove_oid(input_table) #drop oid column is it exists
  check_for_rows(input_table)
  
  # get colnames to process
  if(length(dissolve_list > 0)){
    feature_list <- dissolve_list
  } else{
    feature_list <- colnames(input_table)
  }
  
  saveCount <- 1
  for(col_id in feature_list){
    
    # separate names if combination
    col_ids <- sep_network_names(col_id)
    
    # get list of catchments
    catchments_list <- get_catch_list(col_ids, input_table)
    
    # only proceed if there are catchments in the list. if no catchments, the network will not be included in the output table
    if(length(catchments_list) > 0){
      
      # drop catchments if requested
      if(!is.null(drop_table)){
        drop_list <- get_catch_list(col_ids, drop_table)
        
        catchments_list <- catchments_list[!catchments_list %in% drop_list]
      }
      
      # dissolve based on parameters
      if(calc_area){
        dslv <- catchments_sf %>%
          dplyr::filter(CATCHNUM %in% catchments_list) %>%
          dplyr::summarise(geometry = sf::st_union(geometry)) %>%
          sf::st_buffer(20) %>%
          sf::st_buffer(-20) %>%
          dplyr::mutate(id = col_id,
                        area_km2 = round(as.numeric(sf::st_area(geometry) / 1000000), 2))
      } else{
        dslv <- catchments_sf %>%
          dplyr::filter(CATCHNUM %in% catchments_list) %>%
          dplyr::summarise(geometry = sf::st_union(geometry)) %>%
          sf::st_buffer(20) %>%
          sf::st_buffer(-20) %>%
          dplyr::mutate(id = col_id)
      }
      
      # join AWI if requested
      if(!is.null(intactness_id)){
        if(intactness_id %in% colnames(catchments_sf)){
          awi <- catchments_sf %>%
            dplyr::filter(CATCHNUM %in% catchments_list) %>%
            dplyr::mutate(area = as.numeric(sf::st_area(geometry))) %>%
            sf::st_drop_geometry() %>%
            dplyr::summarise(AWI = sum(.[[intactness_id]] * area) / sum(area))
          
          dslv$AWI <- round(awi$AWI, 4)
        } else{
          warning(paste0("Area-weighted intactness cannot be calculated, ", '"', intactness_id, '"', " not in catchments_sf"))
        }
      }
      
      # set out name
      if(!is.null(out_feature_id)){
        names(dslv)[names(dslv) == "id"] <- out_feature_id
      }
      
      # reorder columns - move geometry to last
      dslv <- dslv  %>%
        dplyr::relocate(geometry, .after = dplyr::last_col())
      
      # append to df
      if(saveCount == 1){
        out_sf <- dslv
        saveCount <- saveCount + 1
      } else{
        out_sf <- rbind(out_sf, dslv)
      }
    }
  }
  return(out_sf)
}
