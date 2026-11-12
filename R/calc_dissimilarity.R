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
#' @param conservation_areas_sf sf object of the features (e.g. conservation areas or networks) to evaluate.
#' @param reserves_id String matching the unique identifier column in \code{conservation_areas_sf}.
#' @param reference_sf sf object of the reference area to compare against.
#' @param raster_layer Single-layer terra SpatRaster that will be clipped to the reference and reserve areas, with crs matching
#' \code{conservation_areas_sf} and \code{reference_sf}.
#' @param raster_type 'continuous' (default) will use the KS-statistic, 'categorical' will use Bray-Curtis.
#' @param categorical_class_values Vector of raster values in \code{raster_layer} to include in the calculation. Required when
#' \code{raster_type = 'categorical'}; values not listed are ignored. Not used for continuous rasters.
#' @param plot_out_dir Path to folder in which to save plots. Default (\code{NULL}) is not to create plots.
#' The folder will be created if it doesn't exist.
#' @param categorical_class_labels Optional data.frame whose first column holds raster values and second column the label to use
#' in Bray-Curtis graphs for each value. Defaults to using the raster values. Labels can be provided for all or a subset of values.
#' @param progress Optional callback function for reporting progress. The function
#'   should accept two arguments: `value`, a numeric value between 0 and 1
#'   indicating the progress, and `detail`, a character string describing the
#'   current processing step. Defaults to `NULL`.
#' @param reference_sample_n Continuous rasters only: maximum number of cells randomly sampled in the reference area.
#' Defaults to 1,000,000.
#' @param reserve_sample_n Continuous rasters only: maximum number of cells randomly sampled in each reserve.
#' Defaults to 500,000.
#' @param seed Random seed set before sampling, so continuous results are reproducible. Defaults to 20260903.
#'
#' @return A numeric vector of dissimilarity values, named by \code{reserves_id} and matching the order of \code{conservation_areas_sf}.
#' Reserves with no valid raster values get \code{NA}. Optionally a dissimilarity plot saved in \code{plot_out_dir} for each
#' computed value.
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
calc_dissimilarity <- function(
    conservation_areas_sf,
    reserves_id,
    reference_sf,
    raster_layer,
    raster_type = c("continuous", "categorical"),
    categorical_class_values = NULL,
    plot_out_dir = NULL,
    categorical_class_labels = data.frame(),
    progress = NULL,
    reference_sample_n = 1000000L,
    reserve_sample_n = 500000L,
    seed = 20260903
) {
  raster_type <- match.arg(raster_type)
  
  report <- function(value, detail) {
    if (!is.null(progress)) {
      progress(value, detail)
    }
  }
  
  # --------------------------------------------------
  # Validate inputs
  # --------------------------------------------------
  reserves_vect <- terra::vect(conservation_areas_sf)
  reference_vect <- terra::vect(reference_sf)
  
  if (!terra::same.crs(reserves_vect, raster_layer)) {
    stop("conservation_areas_sf and raster_layer do not have the same CRS.")
  }
  
  if (!terra::same.crs(reference_vect, raster_layer)) {
    stop("reference_sf and raster_layer do not have the same CRS.")
  }
  
  if (terra::nlyr(raster_layer) != 1L) {
    stop("raster_layer must have exactly one layer.")
  }
  
  if (raster_type == "categorical" &&
      (is.null(categorical_class_values) ||
       length(categorical_class_values) == 0L)) {
    stop(
      "For categorical rasters, provide categorical_class_values."
    )
  }
  
  if (!reserves_id %in% names(conservation_areas_sf)) {
    stop("reserves_id is not a column in conservation_areas_sf.")
  }
  
  # --------------------------------------------------
  # Plot setup
  # --------------------------------------------------
  make_plots <- !is.null(plot_out_dir)
  
  if (make_plots) {
    dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # --------------------------------------------------
  # Helpers
  # --------------------------------------------------
  
  # For continuous rasters: retain a bounded sample, never every cell.
  sample_values_in_polygon <- function(r, v, n) {
    r_crop <- terra::crop(r, terra::ext(v), snap = "out")
    r_mask <- terra::mask(r_crop, v, touches = FALSE)
    
    vals <- terra::spatSample(r_mask, size = n, method = "random", na.rm = TRUE, as.df = TRUE, warn = FALSE)[[1]]
    vals[!is.na(vals)]
  }
  
  # For categorical rasters: retain only one count per class.
  class_counts_in_polygon <- function(r, v, classes) {
    r_crop <- terra::crop(r, terra::ext(v), snap = "out")
    r_mask <- terra::mask(r_crop, v, touches = FALSE)
    
    tab <- terra::freq(r_mask, bylayer = FALSE)
    counts <- setNames(numeric(length(classes)), as.character(classes))
    
    if (!is.null(tab) && nrow(tab) > 0L) {
      keep <- tab$value %in% classes
      
      if (any(keep)) {
        counts[as.character(tab$value[keep])] <- tab$count[keep]
      }
    }
    counts
  }
  
  # Bray-Curtis dissimilarity of class proportions.
  # Normalizing is important: otherwise larger reserves look more
  # dissimilar merely because they contain more raster cells.
  bray_curtis_from_counts <- function(reference_counts, target_counts) {
    if (sum(reference_counts) == 0 || sum(target_counts) == 0) {
      return(NA_real_)
    }
    
    reference_prop <- reference_counts / sum(reference_counts)
    target_prop <- target_counts / sum(target_counts)
    
    sum(abs(reference_prop - target_prop)) / 2
  }
  
  categorical_composition_plot <- function(
    reference_counts,
    target_counts,
    title,
    labels = data.frame()
  ) {
    reference_prop <- reference_counts / sum(reference_counts)
    target_prop <- target_counts / sum(target_counts)
    
    class_names <- names(reference_counts)
    display_names <- class_names
    
    # If labels has at least two columns, treat first as class value
    # and second as display label.
    if (is.data.frame(labels) && ncol(labels) >= 2L) {
      label_index <- match(class_names, as.character(labels[[1]]))
      has_label <- !is.na(label_index)
      
      display_names[has_label] <- as.character(
        labels[[2]][label_index[has_label]]
      )
    }
    
    plot_data <- data.frame(
      class = rep(display_names, 2L),
      proportion = c(reference_prop, target_prop),
      area = rep(
        c("Reference", "Reserve"),
        each = length(display_names)
      )
    )
    
    ggplot2::ggplot(plot_data, ggplot2::aes(x = class, y = proportion, fill = area)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(title = title, x = NULL, y = "Class proportion",  fill = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  # --------------------------------------------------
  # Reference representation, calculated once
  # --------------------------------------------------
  report(0.10, "Preparing reference distribution")
  
  if (raster_type == "continuous") {
    set.seed(seed)
    
    message("Sampling reference values...")
    report(0.20, "Sampling reference values")
    
    reference_vals <- sample_values_in_polygon(raster_layer, reference_vect, n = reference_sample_n)
    
    if (length(reference_vals) == 0L) {
      stop("No non-missing raster values were found in the reference area.")
    }
    
  } else {
    message("Counting reference classes...")
    report(0.20, "Counting reference classes")
    
    reference_counts <- class_counts_in_polygon(raster_layer, reference_vect, classes = categorical_class_values)
    
    if (sum(reference_counts) == 0L) {
      stop("No values matching categorical_class_values were found in the reference area.")
    }
  }
  
  # --------------------------------------------------
  # Main reserve loop
  # --------------------------------------------------
  reserve_ids <- as.character(conservation_areas_sf[[reserves_id]])
  
  result_vector <- numeric(length(reserve_ids))
  names(result_vector) <- reserve_ids
  
  n_reserves <- length(reserve_ids)
  
  if (n_reserves == 0L) {
    return(result_vector)
  }
  
  for (i in seq_along(reserve_ids)) {
    id <- reserve_ids[i]
    
    message(sprintf("[%d/%d] %s", i, n_reserves, id))
    
    report(0.25 + 0.70 * (i - 1) / n_reserves, sprintf("Calculating %d of %d: %s", i, n_reserves, id))
    
    reserve_poly <- reserves_vect[
      reserves_vect[[reserves_id]] == id,
    ]
    
    if (raster_type == "continuous") {
      target_vals <- try(
        sample_values_in_polygon(raster_layer, reserve_poly, n = reserve_sample_n),
        silent = TRUE
      )
      
      if (inherits(target_vals, "try-error") ||
          length(target_vals) == 0L) {
        result_vector[i] <- NA_real_
        next
      }
      
      metric <- ks_stat(reference_vals, target_vals)
      
      if (make_plots) {
        outfile <- file.path(plot_out_dir, paste0(id, ".png"))
        
        plt <- ks_plot(reference_vals, target_vals, plotTitle = paste0(id, " (KS=", round(metric, 3), ")"))
        ggplot2::ggsave(outfile, plt, width = 7, height = 5)
      }
      
    } else {
      target_counts <- try(
        class_counts_in_polygon(raster_layer, reserve_poly, classes = categorical_class_values),
        silent = TRUE
      )
      
      if (inherits(target_counts, "try-error") ||
          sum(target_counts) == 0L) {
        result_vector[i] <- NA_real_
        next
      }
      
      metric <- bray_curtis_from_counts(reference_counts, target_counts)
      
      if (make_plots) {
        outfile <- file.path(plot_out_dir, paste0(id, ".png"))
        
        plt <- categorical_composition_plot(reference_counts, target_counts, title = paste0(id," (Bray-Curtis=", round(metric, 3),")"), labels = categorical_class_labels)
        ggplot2::ggsave(outfile, plt, width = 7, height = 5)
      }
    }
    
    result_vector[i] <- metric
    gc(verbose = FALSE)
  }
  report(1, "Dissimilarity calculation complete")
  result_vector
}