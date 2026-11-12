check_evaluation_table <- function(evaluation_table){
  if(!all(names(evaluation_table) %in% c("class_value", "area_km2", "network", "class_proportion", "target_km2", "prop_target_met"))){
    stop("network_evaluation_table must contain colnames: 'class_value', 'area_km2', 'network', 'class_proportion', 'target_km2', 'prop_target_met' \n Are you using the output from evaluate_targets_using_catchments(), evaluate_targets_using_clip() or evaluate_targets_using_networks()?")
  }
}
