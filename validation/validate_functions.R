# Validation of linear and areal disturbance attributes
# PV 2024-02-02

library(sf)
library(dplyr)
library(summarytools)

validate_project <- function(data_pkg, lyr_line, lyr_poly, output_dir) {

  types <- readr::read_csv('https://raw.githubusercontent.com/beaconsproject/beacons_tools/main/validation/yg_industry_disturbance_types.csv')
  errors <- readr::read_csv('https://raw.githubusercontent.com/beaconsproject/beacons_tools/main/validation/yg_industry_disturbance_types_errors.csv')

  # Create output folder if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Read file and get disturbance attributes
  gpkg <- paste0(data_pkg)
  x_line <- st_read(gpkg, lyr_line, quiet=TRUE) |>
    st_drop_geometry()
  line_vars <- names(x_line)
  x_poly <- st_read(gpkg, lyr_poly, quiet=TRUE) |>
    st_drop_geometry()
  poly_vars <- names(x_poly)

  # Summarize linear attributes
  sink(paste0(output_dir,'/',lyr_line,'_attributes_summary.txt'))
  cat('Project: ', data_pkg, '\n', sep="")
  cat('Date: ', format(Sys.time(), "%d %B %Y"),'\n', sep="")
  cat('\n\n# LINEAR DISTURBANCES\n', sep="")
  for (i in line_vars) {
    cat('\n## Attribute: ',toupper(i),'\n\n', sep="")
    print(dfSummary(x_line[i], graph.col=FALSE, max.distinct.values=20))
  }
  sink()
 
  # Summarize areal attributes
  sink(paste0(output_dir,'/',lyr_poly,'_attributes_summary.txt'))
  cat('Project: ', data_pkg, '\n', sep="")
  cat('Date: ', format(Sys.time(), "%d %B %Y"),'\n', sep="")
  cat('\n\n# AREAL DISTURBANCES\n', sep="")
  for (i in poly_vars) {
    cat('\n## Attribute: ',toupper(i),'\n\n', sep="")
    print(dfSummary(x_poly[i], graph.col=FALSE, max.distinct.values=20))
  }
  sink()

  # Validate linear TYPE_INDUSTRY and TYPE_DISTURBANCE
  line_indu <- types |>
    filter(TYPE_FEATURE=='Linear') |>
    select(TYPE_INDUSTRY) |>
    unique() |>
    pull()
  line_dist <- types |>
    filter(TYPE_FEATURE=='Linear') |>
    select(TYPE_DISTURBANCE) |>
    unique() |>
    pull()
  # These are combinations that are theoretically not permitted but left as is for now
  line_combo <- types |>
    filter(TYPE_FEATURE=='Linear') |>
    select(TYPE_INDUSTRY, TYPE_DISTURBANCE) |>
    mutate(TYPE_COMBINED=paste0(TYPE_INDUSTRY,"***",TYPE_DISTURBANCE)) |>
    unique() |>
    pull()
  # These ones need to be fixed asap - they cannot occur (see fixit script)
  line_combo_error <- errors |>
    filter(TYPE_FEATURE=='Linear') |>
    select(TYPE_INDUSTRY_error, TYPE_DISTURBANCE_error) |>
    mutate(TYPE_COMBINED_error=paste0(TYPE_INDUSTRY_error,"***",TYPE_DISTURBANCE_error)) |>
    unique() |>
    pull()
  x_line_test = mutate(x_line, 
    type_industry_test=ifelse(TYPE_INDUSTRY %in% line_indu, 'ok', 'please fix'),
    type_disturbance_test=ifelse(TYPE_DISTURBANCE %in% line_dist, 'ok', 'please fix'),
    type_combination_test=ifelse(paste0(TYPE_INDUSTRY,"***",TYPE_DISTURBANCE) %in% line_combo, 'ok', 'not expected'),
    type_combination_error=ifelse(!paste0(TYPE_INDUSTRY,"***",TYPE_DISTURBANCE) %in% line_combo_error, 'ok', 'please fix')) |>
    select(TYPE_INDUSTRY, TYPE_DISTURBANCE, type_industry_test, type_disturbance_test, type_combination_test, type_combination_error) |>
    filter(type_industry_test=='please fix' | type_disturbance_test=='please fix' | type_combination_test=='not expected' | type_combination_error=='please fix')
  readr::write_csv(x_line_test, paste0(output_dir, '/',lyr_line,'_types_validation.csv'))

  # Validate areal TYPE_INDUSTRY and TYPE_DISTURBANCE
  poly_indu <- types |>
    filter(TYPE_FEATURE=='Areal') |>
    select(TYPE_INDUSTRY) |>
    unique() |>
    pull()
  poly_dist <- types |>
    filter(TYPE_FEATURE=='Areal') |>
    select(TYPE_DISTURBANCE) |>
    unique() |>
    pull()
  # These are combinations that are theoretically not permitted but left as is for now
  poly_combo <- types |>
    filter(TYPE_FEATURE=='Areal') |>
    select(TYPE_INDUSTRY, TYPE_DISTURBANCE) |>
    mutate(TYPE_COMBINED=paste0(TYPE_INDUSTRY,"***",TYPE_DISTURBANCE)) |>
    unique() |>
    pull()
  x_poly_test = mutate(x_poly, 
    type_industry_test=ifelse(TYPE_INDUSTRY %in% poly_indu, 'ok', 'please fix'),
    type_disturbance_test=ifelse(TYPE_DISTURBANCE %in% poly_dist, 'ok', 'please fix'),
    type_combination_test=ifelse(paste0(TYPE_INDUSTRY,"***",TYPE_DISTURBANCE) %in% poly_combo, 'ok', 'not expected')) |>
    select(TYPE_INDUSTRY, TYPE_DISTURBANCE, type_industry_test, type_disturbance_test, type_combination_test) |>
    filter(type_industry_test=='please fix' | type_disturbance_test=='please fix' | type_combination_test=='not expected')
  readr::write_csv(x_poly_test, paste0(output_dir, '/',lyr_poly,'_types_validation.csv'))

}
