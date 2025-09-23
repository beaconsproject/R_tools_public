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