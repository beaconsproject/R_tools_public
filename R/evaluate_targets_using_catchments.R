### evaluate_targets_using_catchments ###
#
#' Evaluate targets using conservation areas defined by lists of catchments.
#'
#' For a given conservation area or list of conservation area, sums the areas of target values in each conservation area and evaluates 
#' against the targets from gen_targets().
#' Conservation areas are defined by a list of catchment CATCHNUMs, usually produced by the beaconsbuilder package.
#' Usually used for evaluating sets of individual conservation areas.
#' Could also be used to evaluate a network of multiple conservation areas for a set of targets.
#' If evaluating a network of multiple overlapping conservation areas, the network is dissolved (i.e. no double counting of catchments).
#' This can be useful for evaluating the combined area of a network against a set of targets.
#'
#' @param catchments_sf sf object of the catchments dataset with target values summed using [criteria_to_catchments()]
#' @param criteria_name String representing the criteria name that will identify the target area columns (should match criteria name in [criteria_to_catchments()])
#' @param builder_table Data frame where columns are conservation area names and rows are catchments making up the conservation area. 
#' e.g. the "COLUMN_All_Unique_BAs" table output by the beaconsbuilder package.
#' @param target_table Data frame containing columns "class_value", "class_proportion" and "target_km2". i.e. the output from [gen_targets()]. All classes in the target table are 
#'   evaluated. All target classes must have a matching column in the catchments dataset (e.g. class_value 1 in the led target table must have column led_1 in
#'   catchments_sf).
#' @param network_list Vector of networks to evaluate - usually just the list of conservation areas in the builder_table (e.g. PB_0001), but can also be combinations of 
#'   conservation areas (e.g. PB_0001__PB_0002) made using [gen_network_names()]. Usually combinations of conservation areas are multiple overlapping conservation areas that should have 
#'   targets evaluated for their combined area. Combinations of conservation area names must be separated by '__' and the individual respective names must appear in 
#'   builder_table. Defaults to all column names in \code{builder_table}.
#'
#' @return A tibble with columns: 
#'\itemize{
#'  \item{\code{class_value}: the list of class_values from the target_table}
#'  \item{\code{area_km2}: area of the class value in the network}
#'  \item{\code{network}: name of the network copied from network_list}
#'  \item{\code{class_proportion}: copied from target_table}
#'  \item{\code{target_km2}: copied from target_table}
#'  \item{\code{prop_target_met}: area_km2/target_km2}
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' target_table <- gen_targets(ref_poly, led_sample, 1600)
#" # evaluate targets in three individual conservation areas
#' evaluate_targets_using_catchments(
#'   catchments_sample, 
#'   "led", 
#'   builder_table_sample, 
#'   target_table, 
#'   c("PB_0001", "PB_0002", "PB_0003"))
#' 
#' # evaluate targets across a network of multiple conservation areas
#' evaluate_targets_using_catchments(
#'   catchments_sample, 
#'   "led", 
#'   builder_table_sample, 
#'   target_table, 
#'   c("PB_0002__PB_0001"))
evaluate_targets_using_catchments <- function(catchments_sf, criteria_name, builder_table, target_table, network_list=c()){
  
  check_colnames(catchments_sf, cols = "CATCHNUM") # check for CATCHNUM
  check_catchnum_class(catchments_sf, builder_table) # Check catchments match, warning if not
  builder_table <- remove_oid(builder_table) # drop OID column if it exists
  
  # Check all required classes are in the catchments. If not, provide warning. Missing classes will be assumed to have area of zero.
  expected_classes <- paste0(criteria_name, "_", target_table$class_value)
  observed_classes <- as.data.frame(catchments_sf) %>%
    dplyr::select(tidyr::starts_with(criteria_name)) %>%
    names()
  check_classes_in_catchments(expected_classes, observed_classes, warning_ = TRUE, stop_ = FALSE)
  
  # set network_list if not provided
  if(length(network_list) == 0){
    network_list <- colnames(builder_table)
  }
  
  # STEP 1 - sum area of each class in each network
  counter <- 1
  for(net in network_list){
    
    # get catchment list
    nets <- sep_network_names(net)
    
    net_catchments <- get_catch_list(nets, builder_table)
    
    # subset catchments and sum areas
    rows <- catchments_sf %>%
      dplyr::as_tibble() %>%
      dplyr::filter(CATCHNUM %in% net_catchments) %>%
      dplyr::select(CATCHNUM, tidyr::starts_with(criteria_name)) %>%
      tidyr::pivot_longer(cols = tidyr::starts_with(criteria_name), 
                          names_to = "class_value", 
                          names_prefix = paste0(criteria_name, "_"), 
                          values_to = "area", 
                          names_transform = list(class_value = as.integer)) %>%
      dplyr::filter(class_value %in% target_table$class_value) %>% # only sum classes we are evaluating targets for
      dplyr::group_by(class_value) %>%
      dplyr::summarize(area_km2 = sum(area))
    
    missing_classes <- target_table$class_value[!target_table$class_value %in% rows$class_value]
    
    rows <- rows %>%
      dplyr::add_row(class_value = missing_classes, area_km2 = 0) %>%
      dplyr::mutate(network = net) %>%
      dplyr::arrange(class_value)
    
    if(counter == 1){
      df <- rows
      counter <- counter + 1
    } else{
      df <- rbind(df, rows)
    }
  }
  
  # STEP 2 - join targets and calculate proportion of target met
  df <- df %>%
    dplyr::left_join(target_table[c("class_value","class_proportion","target_km2")], by = "class_value")
  
  # add proportion of target met
  df$prop_target_met <- round(df$area_km2 / df$target_km2, 2)
  
  return(df)
}
