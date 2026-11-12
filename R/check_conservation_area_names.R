# check conservation area names are in conservation_area_sf$network
check_conservation_area_names <- function(nets, conservation_areas_sf){
  if(!all(nets %in% conservation_areas_sf$network)){
    stop(paste0("Requested names are not in sf objects 'network' column: ", paste0(nets[!nets %in% conservation_areas_sf$network], collapse=", ")))
  }
}
