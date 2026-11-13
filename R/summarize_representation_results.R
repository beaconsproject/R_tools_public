### summarize_representation_results ###
#
#' Summarize representation results by network.
#'
#' Takes a table of network evaluation results from [evaluate_targets_using_catchments()] or [evaluate_targets_using_clip()]
#' and summarizes either the number of classes passed, or the number of classes not passed (i.e. target gaps) by each network.
#' 
#' Targets of zero are not counted.
#' 
#' Rare targets are often difficult to meet, so \code{target_inclusion_proportion} can be used to drop targets that make up a small 
#' proportion of the total target area. For example \code{target_inclusion_proportion = 0.05} would only consider targets covering 
#' at least 5% of the total target area.
#' 
#' @param network_evaluation_table Data frame output from [evaluate_targets_using_catchments()], [evaluate_targets_using_clip()] or 
#'   [evaluate_targets_using_networks()].
#' @param criteria_name String representing the representation raster. Used to label the summary columns. Usually matches \code{criteria_name} 
#'  when using [evaluate_targets_using_catchments()].
#' @param target_pass_proportion Numeric between 0 and 1 that sets the proportion of the target that needs to be met in order for the target 
#'   to be considered 'passed'. Defaults to 1 (i.e. 100% of the target needs to be met in the network).
#' @param target_inclusion_proportion Numeric between 0 and 1. Target classes with a \code{class_proportion} value less than \code{target_inclusion_proportion} 
#'   are dropped. Defaults to 0 which includes all targets.
#' @param suffix Optional suffix string to add to the end of the summary column name.
#' @param gaps If TRUE, returns a summary column named e.g. \code{led_gaps} with the number of missed targets. Otherwise returns 
#'   a summary column named e.g. \code{led_passed} with the number of passed targets. Defaults to TRUE.
#'
#' @return A tibble with columns \code{network} and a summary column (e.g. \code{led_gaps}).
#' @export
#'
#' @examples
#' target_table <- gen_targets(ref_poly, led_sample, 1600)
#' reserves <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network", 
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' network_evaluation_table <- evaluate_targets_using_clip(
#'   reserves, "network", led_sample, target_table)
#'
#' summarize_representation_results(network_evaluation_table, "led")
#' summarize_representation_results(
#'   network_evaluation_table, "led", 
#'   0.9, 0.05, "_90pcnt_commontargets")
#' summarize_representation_results(network_evaluation_table, "led", 
#'   0.9, suffix = "_90pcnt", gaps = FALSE)
#' 
#' # join multiple summary columns using a list
#' library(purrr)
#' library(dplyr)
#' my_list <- list()
#' my_list <- append(my_list, 
#'   list(summarize_representation_results(network_evaluation_table, "led")))
#' my_list <- append(my_list, 
#'   list(summarize_representation_results(network_evaluation_table, "led", gaps = FALSE)))
#' my_list %>% reduce(left_join, by = "network")
summarize_representation_results <- function(network_evaluation_table, criteria_name, target_pass_proportion = 1, target_inclusion_proportion = 0, suffix = "", gaps = TRUE){
  
  # run columns checks
  missing_cols <- check_colnames(network_evaluation_table, "network")
  if (length(missing_cols) > 0) stop("Column(s) missing from network_evaluation_table: ", paste(missing_cols, collapse = ", "))
  check_evaluation_table(network_evaluation_table)
  
  # make output table template
  df <- dplyr::tibble(network = unique(network_evaluation_table$network))
  
  for(net in df$network){
    
    # drop rows where target is 0
    rslts <- network_evaluation_table[network_evaluation_table$network == net & network_evaluation_table$target_km2 > 0,]
    
    # subset and summarize
    if(gaps){
      # gaps
      df[[paste0(criteria_name,"_gaps", suffix)]][df$network == net] <- nrow(rslts[rslts$prop_target_met < target_pass_proportion & rslts$class_proportion >= target_inclusion_proportion,]) 
    } else{
      #pass count
      df[[paste0(criteria_name,"_passed", suffix)]][df$network == net] <- nrow(rslts[rslts$prop_target_met >= target_pass_proportion & rslts$class_proportion >= target_inclusion_proportion,]) 
    }
  }
  return(df) 
}
