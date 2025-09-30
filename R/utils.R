#check_catchnum <- function(catchments_sf){
  
  # check CATCHNUM exists
#  if(!"CATCHNUM" %in% names(catchments_sf)){
#    stop("Catchments must contain column 'CATCHNUM'")
#  }
#}

check_catchnum_class <- function(catchments_sf, builder_table){
  col_classes <- sapply(colnames(builder_table), function(x) class(builder_table[[x]]))
  if(!all(col_classes == class(catchments_sf$CATCHNUM))){
    warning(paste0("Table column classes do not match class(catchments_sf$CATCHNUM) for columns: ", paste0(colnames(builder_table)[col_classes != class(catchments_sf$CATCHNUM)], collapse=", ")))
  }
}

check_for_geometry <- function(in_sf){
  
  if(!"geometry" %in% names(in_sf)){
    stop("Must contain column: geometry")
  }
}

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

logical_to_integer <- function(x){
  
  if(!x %in% c(TRUE, FALSE)){
    stop("input must be TRUE or FALSE")
  }
  
  return(as.integer(x))
}

check_seeds_areatargets <- function(seeds){
  if(!all(seeds$Areatarget > 0)){
    stop("All Areatarget values must be > 0")
  }
}

check_seeds_in_catchments <- function(seeds, catchments_sf) {
  # Make sure CATCHNUM is character before splitting
  seeds_clean <- seeds %>%
    mutate(CATCHNUM = as.character(CATCHNUM)) %>%
    separate_rows(CATCHNUM, sep = ",") %>%
    mutate(CATCHNUM = trimws(CATCHNUM))  # Remove any whitespace
  
  # Convert to same type as in catchments_sf
  catchnum_vals <- catchments_sf$CATCHNUM
  if (is.numeric(catchnum_vals)) {
    seeds_clean$CATCHNUM <- as.numeric(seeds_clean$CATCHNUM)
  }
  
  # Perform check
  if (!all(seeds_clean$CATCHNUM %in% catchnum_vals)) {
    if (any(seeds_clean$CATCHNUM %in% catchnum_vals)) {
      warning("Not all seeds are in catchments_sf")  # some overlap
    } else {
      stop("None of the seeds are in catchments_sf")  # no overlap
    }
  }
}

check_colnames <- function(x, x_name, cols){
  for(col in cols){
    if(!col %in% colnames(x)){
      stop(paste0("Column '", col, "' not in table '", x_name, "'"))
    }
  }
}

make_all_integer <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.integer(x[[col]])){
        warning(paste0("is.integer(", col, ") == FALSE; converting to integer"))
        x[[col]] <- as.integer(x[[col]])
      }
    }
  }
  return(x)
}

make_all_numeric <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.numeric(x[[col]])){
        warning(paste0("is.numeric(", col, ") == FALSE; converting to numeric"))
        x[[col]] <- as.numeric(x[[col]])
      }
    }
  }
  return(x)
}

make_all_character <- function(x, cols = NULL){
  if(is.null(cols)){
    colss <- colnames(x)
  }else{
    colss <- cols
  }
  for(col in colss){
    if(col %in% colnames(x)){
      if(!is.character(x[[col]])){
        warning(paste0("is.character(", col, ") == FALSE; converting to character"))
        x[[col]] <- as.character(x[[col]])
      }
    }
  }
  return(x)
}

# remove OID from tables comgin out of BUILDER
remove_oid <- function(in_table){
  if("OID" %in% colnames(in_table)){
    out_table <- in_table %>%
      dplyr::select(-OID)
  } else{
    out_table <- in_table
  }
  return(out_table)
}

# Check for rows
check_for_rows <- function(in_table){
  if(nrow(in_table) == 0){
    stop("input_table has no data")
  }
}

# sep_network_names
#' Split network names into their constituent parts.
#'
#' Takes a vector of network names in the form output by [gen_network_names()] (i.e. \code{"__"} separated),
#' and splits the network names back into their original names.
#' 
#' If a single name is provided, a vector of names is returned. If a vector with length >1 is provided, a list
#' of vectors is returned with each list element named with the network name.
#'
#' @param network_names String of a network name, or a vector of multiple network names. Must include the \code{"__"} separator.
#'
#' @return A vector of names, or a list of vectors.
#' @export
#'
#' @examples
#' sep_network_names("PB_1__PB_2")
#' sep_network_names(c("PB_1__PB_2", "PB_1__PB_3"))
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2)
#' sep_network_names(x)
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 4)
#' sep_network_names(x)

sep_network_names <- function(network_names){
  
  if(length(network_names) == 1){
    out_val <- strsplit(network_names, "__")[[1]]
  }
  
  if(length(network_names) > 1){
    out_val <- lapply(network_names, function(x){
      strsplit(x, "__")[[1]]
    })
    names(out_val) <- network_names
  }
  return(out_val)
}

### get_catch_list ###
#
#' Get unique sets of catchment values from a table.
#' 
#' Designed for use with wide format tables such as those output by \code{get_upstream_catchments} where each 
#' column name is a unique protected area id, and each row is a catchment id. Catchment values from all columns 
#' included in CAs_ids are combined and returned after removing duplicates.
#' 
#' @param pa_ids A single column name in input_table, or a vector of multiple column names.
#' @param input_table Data frame where column names are protected area names and rows are catchments making 
#' up each polygon. e.g. output from [get_upstream_catchments()].
#'   
#' @return A vector of unique catchment values.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' get_catch_list("PB_0001", builder_table_sample)
#' get_catch_list(c("PB_0001", "PB_0003"), builder_table_sample)
get_catch_list <- function(CAs_ids, input_table){
  
  # check nets are in input_table
  if(!all(CAs_ids %in% colnames(input_table))){
    stop(paste0("names do not appear as colnames in the input table: ", CAs_ids[!CAs_ids %in% colnames(input_table)]))
  }
  
  # get unique catchments vector from input_table using vector of colnames
  net_catchments <- unique(unlist(lapply(CAs_ids, function(x){
    input_table[[x]]
  })))
  
  net_catchments[!is.na(net_catchments)]
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

### long_to_wide ###
#' Convert long table into wide table with NAs (i.e. BUILDER style table)
#' 
#' Used to prepare long tables for use in BUILDER or functions that take wide tables as 
#' input (e.g.\code{dissolve_catchments_from_table()}).
#' 
#' @param long_df A long style table with values in one column and group names in another.
#' @param col_names The column in long_df from which to get new column names.
#' @param values_col The column in long_df from which to get row values.
#'   
#' @return A wide tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' long_to_wide(
#'   long_df = wide_to_long(builder_table_sample, "group"), 
#'   col_names = "group", 
#'   values_col = "value")
long_to_wide <- function(long_df, col_names, values_col){
  
  # get out table nrow (i.e. longest list of values)
  tbl_rows <- long_df %>%
    dplyr::group_by(.data[[col_names]]) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::summarise(m = max(.data$n)) %>%
    dplyr::pull(.data$m)
  
  values_list <- lapply(unique(long_df[[col_names]]), function(x){
    vals <- long_df[[values_col]][long_df[[col_names]]==x]
    c(vals, rep(NA, tbl_rows - length(vals)))
  })
  names(values_list) <- unique(long_df[[col_names]])
  out_tab <- as.data.frame(do.call(cbind, values_list))
  out_tab <- dplyr::as_tibble(do.call(cbind, values_list))
  
  return(out_tab)
}

# convert list of vectors into df with list element names as colnames and missing values as NAs (i.e. BUILDER style table)
list_to_wide <- function(values_list){
  
  # get out table nrow (i.e. longest list of values)
  tbl_rows <- max(unlist(lapply(values_list, function(x){
    length(x)
  })))
  
  values_list_nas <- lapply(values_list, function(x){
    c(x, rep(NA, tbl_rows - length(x)))
  })
  
  out_tab <- as.data.frame(do.call(cbind, values_list_nas))
  
  return(out_tab)
}

