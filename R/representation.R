### gen_targets ###
#
#' Generate target tables for representation analysis.
#'
#'Takes a reference area polygon and a classified raster for which representation will be evaluated. Calculates the proportion of each class in the reference area 
#'and multiples by reserve size to get target in km2.
#'
#' @param reference_sf sf object of the reference area we are aiming to represent
#' @param representation_raster Raster object of the representation layer classified into categorical classes
#' @param class_values A vector of classes in representation_raster to generate targets for. Defaults to all classes in the representation_raster.
#' @param reserve_size The area in km2 that targets will sum to. Generally the approximate reserve size of reserves being evaluated.
#'
#' @return A tibble with columns: 
#'\itemize{
#'  \item{class_value: the list of class_values}
#'  \item{area_km2: the area of each class_value in the reference_sf}
#'  \item{reserve_size: the provided reserve_size}
#'  \item{class_proportion: area_km2/sum(area_km2)}
#'  \item{target_km2: class_proportion * reserve_area}
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' gen_targets(ref_poly, led_sample, 5000)
gen_targets <- function(reference_sf, representation_raster, reserve_size, class_values = c()){
  
  stopifnot(sf::st_crs(reference_sf) == sf::st_crs(representation_raster))
  
  cellsize <- raster::res(representation_raster)[1] / 1000
  cellarea <- cellsize * cellsize
  
  # if more than one feature in the reference object, dissolve into single geometry
  if(nrow(reference_sf) > 1){
    reference_sf <- reference_sf %>%
      dplyr::summarise(geometry = sf::st_union(.data$geometry))
  }
  
  x <- exactextractr::exact_extract(representation_raster, reference_sf, progress = FALSE)[[1]]
  
  # if no class_values provided, use all values in extracted raster
  if(length(class_values) == 0){
    class_values <- unique(x$value)
    class_values <- class_values[!is.na(class_values)]
  }
  
  # calc sum of each class_value
  df <- x %>%
    dplyr::filter(.data$value %in% class_values) %>%
    dplyr::mutate(area = .data$coverage_fraction * cellarea) %>%
    dplyr::group_by(.data$value) %>%
    dplyr::summarize(area_km2 = sum(.data$area))
  
  names(df) <- c("class_value", "area_km2")
  
  # if any class_values not in df, add them with area of zero
  for(i in class_values){
    if(!i %in% df$class_value){
      df <- rbind(df, data.frame(class_value = i, area_km2 = 0))
    }
  }
  
  # order by class_value
  df <- df[order(df$class_value),]
  
  # add reserve size and calc targets
  df$reserve_size <- reserve_size
  df$class_proportion <- df$area_km2 / sum(df$area_km2)
  df$target_km2 <- round(df$class_proportion * df$reserve_size, 2)
  
  return(df)
}

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


### evaluate_targets_using_clip ###
#
#' Evaluate targets by clipping.
#'
#' Similar to [evaluate_targets_using_catchments()] but instead of pre-calculating the areas of the criteria raster in the catchments,
#' this method clips the criteria raster directly to each provided polygon.
#' Generally used for evaluating conservation areas that did not come from the beaconsbuilder package and do not conform to catchment boundaries.
#' Slower than [evaluate_targets_using_catchments()] for large numbers of polygons.

#' @param conservation_areas_sf sf object containing conservation areas to evaluate
#' @param conservation_areas_id Unique id column in conservation_areas_sf containing conservation area names as strings.
#' @param representation_raster Raster object of the criteria layer that will be evaluated, with crs matching conservation_areas_sf
#' @param target_table Data frame containing columns "class_value", "class_proportion" and "target_km2". i.e. the output from [gen_targets()]. All classes in the target table are 
#'   evaluated. All class_values must match the values in the representation_raster
#' @param conservation_areas_list The conservation_areas_id values to process. Defaults to all conservation_areas_id's in conservation_areas_sf
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
#' pas <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample,
#'   "network",
#'   dissolve_list = c("PB_0001", "PB_0002"))
#' evaluate_targets_using_clip(pas, "network", led_sample, target_table)
evaluate_targets_using_clip <- function(conservation_areas_sf, conservation_areas_id, representation_raster, target_table, conservation_areas_list=c()){
  
  stopifnot(sf::st_crs(conservation_areas_sf) == sf::st_crs(representation_raster))
  
  # set conservation_areas_list if not provided
  if(length(conservation_areas_list) == 0){
    conservation_areas_list <- conservation_areas_sf[[conservation_areas_id]]
  }
  
  cell_area <- prod(raster::res(representation_raster)) / 1000000 # convert to area in km2, assumes raster res is in metres
  
  counter <- 1
  for(reserve in conservation_areas_list){
    
    # STEP 1 - sum area of each class in each network
    # subset
    reserve_sf <- conservation_areas_sf[conservation_areas_sf[[conservation_areas_id]] == reserve,]
    
    # clip raster
    x <- exactextractr::exact_extract(representation_raster, reserve_sf, progress = FALSE)[[1]]
    names(x) <- c("class_value", "coverage_fraction")
    rows <- x %>%
      dplyr::filter(.data$class_value %in% target_table$class_value) %>%
      dplyr::mutate(area = .data$coverage_fraction * cell_area) %>%
      dplyr::group_by(.data$class_value) %>%
      dplyr::summarize(area_km2 = sum(.data$area))
    
    missing_class_values <- target_table$class_value[!target_table$class_value %in% rows$class_value]
    
    rows <- rows %>% 
      dplyr::add_row(class_value = missing_class_values, area_km2 = 0) %>%
      dplyr::mutate(network = reserve) %>%
      dplyr::arrange(.data$class_value)
    
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

### evaluate_criteria_using_clip ###
#
#' Evaluate criteria by clipping with a conservation areas or a network.
#'
#' Clips the criteria raster to the conservation area and calculates the proportion of each class within it. If a unique 
#' identifer (conservation_area_id) is provided, proportions are computed separately for each corresponding polygon. 
#'
#' @param conservation_area_sf sf object of the reference area we are aiming to represent
#' @param conservation_area_id Unique id column in conservation_areas_sf containing conservation area names as strings. Can be a unique identifier representing a network. 
#' @param criteria_raster Raster object of the representation layer classified into categorical classes
#' @param class_values A vector of classes in representation_raster to generate targets for. Defaults to all classes in the representation_raster.
#' @param target_size The area in km2 that targets will sum to. Generally the approximate target size of reference area being evaluated. set target Must be set to TRUE

#' @return A tibble with columns: 
#'\itemize{
#'  \item{class_value: the list of class_values}
#'  \item{area_km2: the area of each class_value in the conservation_area_sf}
#'  \item{class_proportion: area_km2/sum(area_km2)}
#'  
#'  }
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' evaluate_criteria_using_clip(ref_poly, led_sample)
evaluate_criteria_using_clip <- function(conservation_area_sf, criteria_raster, conservation_area_id = NULL, class_values = c(), target_size = NULL){
  
  stopifnot(sf::st_crs(conservation_area_sf) == sf::st_crs(criteria_raster))

  ## Dissolve geometry
  if (!is.null(conservation_area_id) && conservation_area_id %in% names(conservation_area_sf)) {
    conservation_area_sf <- conservation_area_sf |>
      dplyr::group_by(.data[[conservation_area_id]]) |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
  #} else if (nrow(conservation_area_sf) > 1L) {
  } else {
    conservation_area_sf <- conservation_area_sf |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
  }
  
  cellsize <- raster::res(criteria_raster)[1] / 1000
  cellarea <- cellsize * cellsize
  
  x <- exactextractr::exact_extract(criteria_raster, conservation_area_sf, include_cols = conservation_area_id, progress = FALSE)
  
  x_df <- purrr::map_dfr(x, ~ .)
  
  # if no class_values provided, use all values in extracted raster
  if(length(class_values) == 0){
    class_values <- unique(x_df$value)
    class_values <- class_values[!is.na(class_values)]
  }
  
  group_vars <- if (!is.null(conservation_area_id) && conservation_area_id != "") c(conservation_area_id, "value") else "value"
  
  # calc sum of each class_value
  df <- x_df %>%
    dplyr::filter(value %in% class_values) %>%
    dplyr::mutate(area = coverage_fraction * cellarea) %>%
    dplyr::group_by_at(group_vars) %>%
    #dplyr::group_by(!!rlang::sym(conservation_area_id), value) %>%
    dplyr::summarize(area_km2 = sum(area), .groups = "drop")
  
  names(df) <- c(conservation_area_id, "class_value", "area_km2")
  
  # if any class_values not in df, add them with area of zero
  for(i in class_values){
    if(!i %in% df$class_value){
      df <- rbind(df, data.frame(class_value = i, area_km2 = 0))
    }
  }
  
  # order by class_value and calculate class_proportion
  if (!is.null(conservation_area_id) && conservation_area_id != "") {
    df <- df %>%
      dplyr::arrange(!!rlang::sym(conservation_area_id), class_value) %>%
      dplyr::group_by(!!rlang::sym(conservation_area_id)) %>%
      dplyr::mutate(class_proportion = area_km2 / sum(area_km2)) %>%
      dplyr::ungroup()
  } else {
    df <- df %>%
      dplyr::arrange(class_value) %>%
      dplyr::mutate(class_proportion = area_km2 / sum(area_km2))
  }
  
  if(!is.null(target_size)){
    df$target_size <- target_size
    df$target_km2 <- round(df$class_proportion * df$target_size, 2)
  }
  
  return(df)
}

### calc_dissimilarity ###
#
#' Calculate dissimilarity values between a set of polygons and a reference area.
#'
#' For a list of features (e.g. conservation areas or networks), calculate the dissimilarity value between each feature and a reference area for the 
#' provided raster layer. Continuous rasters use the KS-statistic to compare distributions, categorical rasters use the Bray-Curtis
#' statistic. Graphs comparing distributions can optionally be created and saved to a user provided file path.
#' 
#' NA values are always removed. For categorical rasters, values can optionally be subset for the calculation and graphs using 
#' \code{categorical_class_values}.
#'
#' @param reserves_sf sf object with unique id column named \code{network}
#' @param reference_sf sf object of the reference area to compare against.
#' @param raster_layer Raster object that will be clipped to the reference and reserve areas, with crs matching reserves_sf
#' @param raster_type 'categorical' will use Bray-Curtis, 'continuous' will use KS-statistic.
#' @param categorical_class_values Vector of raster values in \code{raster_layer} of type 'categorical' to include in the calculation. 
#' Allows unwanted values to be dropped. Defaults to include all non-NA values.
#' @param plot_out_dir Path to folder in which to save plots. Default is not to create plots. 
#' Only creates plots if valid file path is provided. Dir will be created if it doesn't exist.
#' @param categorical_class_labels Optional data.frame object with columns \code{values} and \code{labels} indicating the label to use in Bray-Curtis graphs
#' for each raster value. Defaults to using the raster values. Labels can be provided for all or a subset of values. See examples.
#'
#' @return A vector of dissimilarity values matching the order of the input \code{reserves_sf}. Optionally a dissimilarity plot saved in the 
#' \code{plot_out_dir} for each computed value.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats ks.test
#' @export
#'
#' @examples
#' reserves <- dissolve_catchments_from_table(
#'   catchments_sample, 
#'   builder_table_sample, 
#'   "network",
#'   dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
#' calc_dissimilarity(reserves, "network", ref_poly, led_sample, 'categorical')
#' calc_dissimilarity(reserves, , ref_poly, led_sample, 'categorical', c(1,2,3,4,5))
#' calc_dissimilarity(reserves, "network", ref_poly, led_sample, 'categorical', c(1,2,3,4,5), 
#'   "C:/temp/plots", data.frame(values=c(1,2,3,4,5), labels=c("one","two","three","four","five")))
#' calc_dissimilarity(reserves, ref_poly, led_sample, 'continuous', plot_out_dir="C:/temp/plots")
calc_dissimilarity <- function(reserves_sf, reserves_id, reference_sf, raster_layer, raster_type, categorical_class_values=c(), plot_out_dir=NULL, categorical_class_labels=data.frame()){
  
  # geometries should match
  stopifnot(sf::st_crs(reserves_sf) == sf::st_crs(reference_sf))
  stopifnot(sf::st_crs(reserves_sf) == sf::st_crs(raster_layer))
  
  # set up output vector
  result_vector <- c()
  
  # check raster_layer is valid
  if(!raster_type %in% c('categorical', 'continuous')){
    stop("raster_layer must be on of: 'categorical', 'continuous'")
  }
  
  # check geometry column is present in sf objects
  check_for_geometry(reference_sf)
  check_for_geometry(reserves_sf)
  
  # should plots be made? Attempt to create directory if it doesn't already exist
  if(!is.null(plot_out_dir)){
    make_plots <- TRUE
    dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)
  } else{
    make_plots <- FALSE
  }
  
  # get the reference area and make sure it is a single geometry
  ref_sf <- reference_sf %>%
    dplyr::summarise(geometry = sf::st_union(.data$geometry))
  
  # extract values in ref area
  ref_ext <- exactextractr::exact_extract(raster_layer, reference_sf, progress = FALSE)[[1]]
  
  # filter by class_vals if provided, remove NAs, only keep cells with majority in reference_sf
  if(raster_type == 'categorical' & length(categorical_class_values) > 0){
    reference_vals <- ref_ext %>%
      dplyr::filter(!is.na(.data$value)) %>%
      dplyr::filter(.data$value %in% categorical_class_values) %>%
      dplyr::filter(.data$coverage_fraction > 0.5) %>% # only keep values from cells with at least half their area in the polygon
      dplyr::pull(.data$value)
  } else{
    reference_vals <- ref_ext %>%
      dplyr::filter(!is.na(.data$value)) %>%
      dplyr::filter(.data$coverage_fraction > 0.5) %>%
      dplyr::pull(.data$value)
  }
  
  # split reserve into blocks of 10 for processing
  net_list <- dplyr::pull(reserves_sf, !!reserves_id) |> as.character()
  net_list_grouped <- split(net_list, ceiling(seq_along(net_list)/10))
  
  # run in blocks of 10 - seems optimal for maintaining a fast extract
  counter <- 1
  for(net_list_g in net_list_grouped){
    
    message(paste0("processing ", counter*10-9, " of ", length(net_list)))
    counter <- counter + 1
    
    reserves_sf_g <- reserves_sf %>%
      dplyr::filter(!!rlang::sym(reserves_id) %in% net_list_g) # subset dissolved networks by the block of networks
    x <- exactextractr::exact_extract(raster_layer, reserves_sf_g, progress = FALSE) # extract
    
    names(x) <- net_list_g # name the list elements by their associated netname
    
    for(net in net_list_g){
      
      # for each network in the block, extract the values...
      if(raster_type == 'categorical' & length(categorical_class_values) > 0){
        target_vals <- x[[net]] %>%
          dplyr::filter(!is.na(.data$value)) %>%
          dplyr::filter(.data$value %in% categorical_class_values) %>%
          dplyr::filter(.data$coverage_fraction > 0.5) %>% # only keep values from cells with at least half their area in the polygon
          dplyr::pull(.data$value)
      } else{
        target_vals <- x[[net]] %>%
          dplyr::filter(!is.na(.data$value)) %>%
          dplyr::filter(.data$coverage_fraction > 0.5) %>%
          dplyr::pull(.data$value)
      }
      
      # run dissimilarity
      if(raster_type == "categorical"){
        result <- bc_stat(reference_vals, target_vals)
      } else{
        result <- ks_stat(reference_vals, target_vals)
      }
      
      # add result t return vector
      result_vector <- c(result_vector, result)
      
      # generate plot if requested
      if(make_plots){
        plot_out_path <- file.path(plot_out_dir, paste0(net, ".png"))
        
        if(raster_type == 'categorical'){
          plt <- bc_plot(reference_vals, target_vals, plotTitle = paste0(net, " (BC=", result, ")"), labels=categorical_class_labels)
        } else{
          plt <- ks_plot(reference_vals, target_vals, plotTitle = paste0(net, " (KS=", result, ")"))
        }
        ggplot2::ggsave(plt, file=plot_out_path)
      }
    }
  }
  return(result_vector)
}

# This should be written in RCPP
ks_stat <- function(refVal, netVal) {
  # calculate KS statistic (representation index)
  ri <- suppressWarnings(round(ks.test(refVal, netVal)[[1]][[1]], 3))
  return(ri)
}

# This should be written in RCPP
bc_stat <- function(refVal, netVal) {
  
  x1 <- dplyr::as_tibble(refVal) %>%
    dplyr::count(.data$value)
  names(x1) <- c("cat","strata")
  
  x2 <- dplyr::as_tibble(netVal) %>%
    dplyr::count(.data$value)
  names(x2) <- c("cat","reserve")
  
  x <- merge(x1,x2,by="cat",all=T)
  #x$cat <- as.character(x$cat)
  x$strata <- as.numeric(x$strata)
  x$reserve <- as.numeric(x$reserve)
  x$reserve[is.na(x$reserve)] <- 0
  x$reserve[is.na(x$strata)] <- 0 # this is needed in case there is one reserve pixel and no strata pixel
  x$strata[is.na(x$strata)] <- 0
  x$strata <- x$strata/sum(x$strata)
  x$reserve <- x$reserve/sum(x$reserve)
  
  # calculate Bray-Curtis dissimilariy
  ri <- round(sum(abs(x$strata-x$reserve))/(sum(x$strata)+sum(x$reserve)), 3)
  return(ri)
}

ks_plot <- function(refVal, netVal, plotTitle="") {
  
  regLab <- "Reference area"
  netLab <- "Network"
  
  z1 <- c(refVal, netVal)
  z2 <- c(rep(regLab,length(refVal)), rep(netLab,length(netVal)))
  zz <- data.frame(cbind(z1,z2),stringsAsFactors=FALSE)
  names(zz) <- c("values","criteria")
  zz$values <- round(as.numeric(zz$values),3)
  
  # create and save density plot
  p <- ggplot2::ggplot(zz, ggplot2::aes(x=.data$values)) + ggplot2::geom_density(ggplot2::aes(group=.data$criteria, color=.data$criteria)) +
    ggplot2::ggtitle(plotTitle) +
    ggplot2::labs(x="Indicator value", y="Density")
  
  return(p)
}

bc_plot <- function(refVal, netVal, plotTitle="", labels=data.frame()) {
  
  x1 <- dplyr::as_tibble(refVal) %>%
    dplyr::count(.data$value)
  names(x1) <- c("cat","strata")
  
  x2 <- dplyr::as_tibble(netVal) %>%
    dplyr::count(.data$value)
  names(x2) <- c("cat","reserve")
  
  x <- merge(x1,x2,by="cat",all=T)
  x <- x[order(as.integer(as.character(x$cat))),]
  x$strata <- as.numeric(x$strata)
  x$reserve <- as.numeric(x$reserve)
  x$reserve[is.na(x$reserve)] <- 0
  x$reserve[is.na(x$strata)] <- 0 # this is needed in case there is one reserve pixel and no strata pixel
  x$strata[is.na(x$strata)] <- 0
  x$strata <- x$strata/sum(x$strata) #as.integer(x$strata)
  x$reserve <- x$reserve/sum(x$reserve) #as.integer(x$reserve)
  
  # prep labels if present
  if(nrow(labels) > 0 & "values" %in% names(labels)){
    for(i in labels$values){
      if(i %in% x$cat){
        x$cat[x$cat == i] <- labels$label[labels$values == i]
      }
    }
  }
  x$cat <- factor(x$cat, levels = x$cat)
  
  p <- ggplot2::ggplot(x, ggplot2::aes(x=.data$cat, y=.data$reserve)) + ggplot2::geom_bar(stat="identity", fill="white", colour="black") + ggplot2::coord_flip()
  p <- p + ggplot2::geom_point(data=x, ggplot2::aes(x=.data$cat, y=.data$strata), colour="black", size=3) + ggplot2::theme(legend.position = "none")
  p <- p + ggplot2::labs(x="", y="Proportional area (dots indicate regional proportions)")
  p <- p + ggplot2::ggtitle(plotTitle)
  
  return(p)
}

#' Sum areas of criteria raster values inside catchment polygons.
#'
#' For a given raster layer, sums the area of all unique values in the raster and adds summed area's as columns in the catchments dataset.
#' 
#' Raster projection units are assumed to be in metres.
#'
#' @param catchments_sf sf object of catchments
#' @param criteria_raster Raster object of the criteria layer that will be summed, with crs matching catchments
#' @param criteria_name String representing the criteria name that will provide the suffix in the new column names
#' @param class_vals Vector of the class values to sum. Defaults to including all unique values of the raster intersecting the catchments
#'
#' @return sf object matching catchments_sf, with the additional columns added and areas reported in km2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' criteria_to_catchments(catchments_sample, led_sample, "led")

criteria_to_catchments <- function(catchments_sf, criteria_raster, criteria_name, class_vals = c()){
  
  check_colnames(catchments_sf, cols = "CATCHNUM") # check for CATCHNUM
  stopifnot(sf::st_crs(catchments_sf) == sf::st_crs(criteria_raster))
  
  cell_area <- prod(raster::res(criteria_raster)) / 1000000 # convert to area in km2, assumes raster res is in metres
  
  # split catchments into blocks of 50 for processing
  #catch_list <- unique(as.character(catchments_sf$CATCHNUM))
  catch_list <- unique(catchments_sf$CATCHNUM)
  catch_list_grouped <- split(catch_list, ceiling(seq_along(catch_list)/50))
  
  block_counter <- 1
  catch_counter <- 1
  for(catch_list_i in catch_list_grouped){
    
    message(paste0("block ", block_counter, " of ", length(catch_list_grouped)))
    block_counter <- block_counter + 1
    
    catchments_i <- catchments_sf[catchments_sf$CATCHNUM %in% catch_list_i,] # subset catchments
    x <- exactextractr::exact_extract(criteria_raster, catchments_i, progress = FALSE) # extract
    
    names(x) <- catch_list_i # name the list elements by their associated CATCHNUM. Catchnum gets converted to character by names().
    
    # sum all values areas. Filter by class_vals later when we can calculate all unique values in all catchments
    for(catch_i in catch_list_i){
      i_sums <- x[[as.character(catch_i)]] %>%
        dplyr::mutate(area = .data$coverage_fraction * cell_area) %>%
        dplyr::group_by(.data$value) %>%
        dplyr::summarise(area_km2 = sum(.data$area)) %>%
        dplyr::mutate(CATCHNUM = catch_i)
      
      # append into long tibble
      if(catch_counter == 1){
        df_long <- i_sums
        catch_counter <- catch_counter + 1
      } else{
        df_long <- rbind(df_long, i_sums)
      }
    }
  }
  
  # set up class_vals
  df_long_vals <- unique(df_long$value)
  if(length(class_vals) == 0){
    class_vals <- df_long_vals
  }
  
  missing_class_vals <- class_vals[!class_vals %in% df_long_vals]
  missing_values <- as.character(missing_class_vals[[1]])
  
  # pivot to wide table
  df_wide <- df_long %>%
    mutate(
      CATCHNUM = CATCHNUM,
      value = as.character(value),
      area_km2 = as.numeric(area_km2)
    ) %>%
    complete(
      CATCHNUM,
      value = as.character(missing_values),
      fill = list(area_km2 = 0)
    ) %>%
    pivot_wider(
      id_cols = CATCHNUM,
      names_from = value,
      values_from = area_km2,
      names_prefix = paste0(criteria_name, "_")
    )

  catchments_sf <- catchments_sf %>%
    dplyr::select(-dplyr::matches(setdiff(names(df_wide), "CATCHNUM"))) %>% # remove any df_wide columns already in catchments_sf. Effectively overwrites the old with the new columns
    dplyr::left_join(df_wide, by = "CATCHNUM") %>%
    dplyr::select(!dplyr::ends_with("NA"))
  
  return(catchments_sf)
}

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
  check_colnames(network_evaluation_table, "evaluation_table", cols = "network")
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
