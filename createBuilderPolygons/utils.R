prep_input_column <- function(out_dir, type= "BENCHMARKS"){
  if(is.null(out_dir)){
    stop("You need to point to a Builder_output directory")
  }
  
  if(type =="BENCHMARKS"){
    benchmarks_out <- list.files(out_dir, pattern = "COLUMN_All_Unique_BAs.csv")
  }else if (type == "UPSTREAM"){
    benchmarks_out <- list.files(outdir, pattern = "UPSTREAM_CATCHMENTS_COLUMN.csv")
  }else if (type == "DOWNSTREAM"){
    benchmarks_out <- list.files(outdir, pattern = "DOWNSTREAM_CATCHMENTS_COLUMN.csv")
  }else{
    stop("Wrong Type. Please select BENCHMARK, UPSTREAM or DOWNSTREAM")
  }
  
  if(length(benchmarks_out) > 1){
    benchmarks_out <- benchmarks_out[[length(benchmarks_out)]]
    warning("Multiple output tables in output folder, returning newest")
  } else if (length(benchmarks_out)==0){ # If no new table was added, for out_dir, do nothing. For temp dir, delete. Then throw error.
    stop("No output produced. Check parameters, file paths and input tables.")
  }
  
  message(paste0("Returning BUILDER table: ", benchmarks_out))
  out_tab <- utils::read.csv(file.path(outdir, benchmarks_out))
  
  # make sure BUILDER made some conservation areas
  if(nrow(out_tab) == 0){
    warning("No conservation areas to return, try changing the BUILDER parameters by decreasing intactness and/or area requirements")
    return(NULL)
  }
  
  # return
  out_tab$OID <- NULL
  return(dplyr::as_tibble(out_tab))
}
