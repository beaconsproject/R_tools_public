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
