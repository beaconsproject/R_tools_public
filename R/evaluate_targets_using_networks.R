### evaluate_targets_using_networks ###
#
#' Combines individual conservation area results into networks.
#'
#' For a set of conservation areas in a network, tests if at least one conservation area met the target.
#' 
#' Each target therefore has to be met in full by at least one conservation area in the network. To assess target classes summed 
#' across the entire network, use [evaluate_targets_using_catchments()] or [evaluate_targets_using_clip()].
#' 
#' @param individual_results Tibble or data frame output from [evaluate_targets_using_catchments()] or [evaluate_targets_using_clip()]
#' @param network_list Vector of networks to evaluate, such as those made with [gen_network_names()].
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
#' results <- evaluate_targets_using_catchments(
#'   catchments_sample, 
#'   "led", 
#'   builder_table_sample, 
#'   target_table, 
#'   colnames(builder_table_sample))
#'  
#' # evaluate networks of 2 conservation areas 
#' network_list <- gen_network_names(unique(results$network), 2)
#' evaluate_targets_using_networks(results, network_list)
#' 
#' # evaluate the option of using all three conservation areas in a single network
#' evaluate_targets_using_networks(results, "PB_0001__PB_0002__PB_0003")
evaluate_targets_using_networks <- function(individual_results, network_list){
  
  missing_cols <- check_colnames(individual_results, "network")
  if (length(missing_cols) > 0) stop("Column(s) missing from individual_results: ", paste(missing_cols, collapse = ", "))
  individual_results$network <- as.character(individual_results$network)
  check_evaluation_table(individual_results)
  
  # process 1000 networks at a time to limit the number of bin_row calls
  network_list_grouped <- split(network_list, ceiling(seq_along(network_list)/1000))
  counter <- 1
  for(network_group in network_list_grouped){
    
    group_list <- list()
    
    for(net in network_group){
      
      nets <- sep_network_names(net)
      
      check_conservation_area_names(nets, individual_results)
      
      # rebuild results table with whichever benchmark has the largest area for each target
      net_rslts <- individual_results %>%
        dplyr::filter(.data$network %in% nets, .data$target_km2 > 0) %>% # drop targets that are zero
        dplyr::group_by(.data$class_value) %>%
        dplyr::arrange(dplyr::desc(.data$area_km2)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::mutate(network = net)
      group_list <- append(group_list, list(net_rslts))
    }
    group_df <- dplyr::bind_rows(group_list) # bind group into df
    
    if(counter == 1){
      df <- group_df
    } else{
      df <- dplyr::bind_rows(df, group_df)
    }
    counter <- counter + 1
  }
  return(df)
}
