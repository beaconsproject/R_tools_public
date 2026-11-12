# check all required criteria classes have been added to catchments for a given input target table
check_classes_in_catchments <- function(expected_classes, observed_classes, warning_, stop_){
  
  missing_classes <- expected_classes[!expected_classes %in% observed_classes]
  
  if(warning_ & length(missing_classes) > 0){
    warning(paste0("Classes are not in catchments: ", paste0(missing_classes, collapse=", "), ". Add them to catchments using criteria_to_catchments(), otherwise an area an area of zero will be assumed."))
  }
  if(stop_  & length(missing_classes) > 0){
    stop(paste0("Classes are not in catchments: ", paste0(missing_classes, collapse=", "), ". Add them to catchments using criteria_to_catchments()."))
  }
}
