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
    
    reserves_sf_g <- dplyr::filter(reserves_sf, !!reserves_id %in% net_list_g) # subset dissolved networks by the block of networks
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