
# seeds
test_that("seeds uses all catchments with a single area target", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  expect_s3_class(s, "tbl_df")
  expect_named(s, c("CATCHNUM", "Areatarget"))
  expect_equal(nrow(s), nrow(catchments_sample))
  expect_type(s$CATCHNUM, "integer")
  expect_true(all(s$Areatarget == 1e9))
})

test_that("seeds takes area targets from a column", {
  catch <- catchments_sample
  catch$area_target <- 5e8
  s <- seeds(catch, areatarget_col = "area_target")
  expect_equal(nrow(s), nrow(catch))
  expect_true(all(s$Areatarget == 5e8))
})

test_that("seeds single value has priority over column", {
  catch <- catchments_sample
  catch$area_target <- 5e8
  s <- seeds(catch, areatarget_value = 1e9, areatarget_col = "area_target")
  expect_true(all(s$Areatarget == 1e9))
})

test_that("seeds filters by polygon and joins area targets from polygon", {
  # 15 km circle around the centre of the sample catchments
  poly <- st_sf(
    Areatarget = 2e9,
    geometry = st_buffer(st_centroid(st_union(catchments_sample)), 15000)
  )
  s <- seeds(catchments_sample, filter_polygon = poly,
             areatarget_polygon = poly, areatarget_polygon_col = "Areatarget")
  expect_equal(nrow(s), 15)
  expect_true(all(s$CATCHNUM %in% catchments_sample$CATCHNUM))
  expect_true(all(s$Areatarget == 2e9))
})

test_that("seeds errors when no area target is provided", {
  expect_error(seeds(catchments_sample), "No valid area target method provided")
})

test_that("seeds warns and errors when area target column is missing", {
  expect_error(
    expect_warning(seeds(catchments_sample, areatarget_col = "not_a_col"), "not_a_col"),
    "No valid area target method provided"
  )
})

test_that("seeds errors when polygon area targets are not numeric", {
  poly <- st_sf(
    Areatarget = "big",
    geometry = st_buffer(st_centroid(st_union(catchments_sample)), 15000)
  )
  expect_error(
    seeds(catchments_sample, areatarget_polygon = poly, areatarget_polygon_col = "Areatarget"),
    "must be numeric"
  )
})

test_that("seeds errors when filter polygon selects no catchments", {
  # polygon far away from the sample catchments
  poly <- st_sf(geometry = st_buffer(st_sfc(st_point(c(0, 0)), crs = st_crs(catchments_sample)), 1000))
  expect_error(
    seeds(catchments_sample, filter_polygon = poly, areatarget_value = 1e9),
    "No catchments selected"
  )
})


# neighbours
test_that("neighbours returns expected structure", {
  n <- neighbours(catchments_sample)
  expect_s3_class(n, "tbl_df")
  expect_named(n, c("CATCHNUM", "neighbours", "key"))
  expect_equal(nrow(n), 446)
  expect_type(n$CATCHNUM, "integer")
  expect_type(n$neighbours, "integer")
  expect_equal(n$key, 0:(nrow(n) - 1))
})

test_that("neighbours never lists a catchment as its own neighbour", {
  n <- neighbours(catchments_sample)
  expect_false(any(n$CATCHNUM == n$neighbours))
})

test_that("neighbours are symmetric", {
  # if A neighbours B, B must neighbour A
  n <- neighbours(catchments_sample)
  pairs    <- paste(n$CATCHNUM, n$neighbours)
  reversed <- paste(n$neighbours, n$CATCHNUM)
  expect_setequal(pairs, reversed)
})

test_that("neighbours errors on a longitude/latitude CRS", {
  catch_ll <- st_transform(catchments_sample[1:5, ], 4326)
  expect_error(neighbours(catch_ll), "projected CRS")
})

test_that("neighbours only uses catchments from the input", {
  n <- neighbours(catchments_sample)
  expect_true(all(n$CATCHNUM %in% catchments_sample$CATCHNUM))
  expect_true(all(n$neighbours %in% catchments_sample$CATCHNUM))
})


# reserve_seeds
test_that("reserve_seeds returns one row per conservation area", {
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001", "PB_0002"))
  r <- reserve_seeds(catchments_sample, cas, "network", areatarget_value = 1e9)
  expect_s3_class(r, "tbl_df")
  expect_named(r, c("network", "Areatarget", "CATCHNUM"))
  expect_equal(r$network, c("PB_0001", "PB_0002"))
  expect_true(all(r$Areatarget == 1e9))
})

test_that("reserve_seeds returns integer Areatarget, like seeds", {
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001"))
  r <- reserve_seeds(catchments_sample, cas, "network", areatarget_value = 1e9)
  expect_type(r$Areatarget, "integer")
})

test_that("reserve_seeds errors on invalid joinType", {
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001"))
  expect_error(reserve_seeds(catchments_sample, cas, "network", 1e9, joinType = "WITHIN"), "should be one of")
})

test_that("reserve_seeds does not change the scipen option", {
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001"))
  withr::local_options(scipen = 0)
  reserve_seeds(catchments_sample, cas, "network", 1e9)
  expect_equal(getOption("scipen"), 0)
})

test_that("reserve_seeds CENTROID returns the catchments of each conservation area", {
  # conservation areas built from builder_table_sample, so CENTROID join
  # should give back exactly the catchments listed in the table
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001", "PB_0002"))
  r <- reserve_seeds(catchments_sample, cas, "network", areatarget_value = 1e9, joinType = "CENTROID")
  for (i in seq_len(nrow(r))) {
    expect_setequal(
      as.integer(strsplit(r$CATCHNUM[i], ",")[[1]]),
      get_catch_list(r$network[i], builder_table_sample)
    )
  }
})

test_that("reserve_seeds INTERSECT includes all CENTROID catchments", {
  # INTERSECT also picks up catchments touching the edge of the conservation area
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001", "PB_0002"))
  r_cent <- reserve_seeds(catchments_sample, cas, "network", 1e9, joinType = "CENTROID")
  r_int  <- reserve_seeds(catchments_sample, cas, "network", 1e9, joinType = "INTERSECT")
  for (i in seq_len(nrow(r_cent))) {
    cent <- strsplit(r_cent$CATCHNUM[i], ",")[[1]]
    int  <- strsplit(r_int$CATCHNUM[i], ",")[[1]]
    expect_true(all(cent %in% int))
    expect_gte(length(int), length(cent))
  }
})


# check_seeds_areatargets
test_that("check_seeds_areatargets passes when all area targets are > 0", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  expect_silent(check_seeds_areatargets(s))
})

test_that("check_seeds_areatargets errors on zero or negative area targets", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  s$Areatarget[1] <- 0
  expect_error(check_seeds_areatargets(s), "All Areatarget values must be > 0")
  s$Areatarget[1] <- -5
  expect_error(check_seeds_areatargets(s), "All Areatarget values must be > 0")
})


# check_seeds_in_catchments
test_that("check_seeds_in_catchments passes when all seeds are in catchments", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  expect_silent(check_seeds_in_catchments(s, catchments_sample))
})

test_that("check_seeds_in_catchments accepts comma separated CATCHNUMs from reserve_seeds", {
  cas <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network",
                                        dissolve_list = c("PB_0001", "PB_0002"))
  r <- reserve_seeds(catchments_sample, cas, "network", 1e9)
  expect_silent(check_seeds_in_catchments(r, catchments_sample))
})

test_that("check_seeds_in_catchments warns when some seeds are missing", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  s$CATCHNUM[1] <- 1L # not a CATCHNUM in the sample
  expect_warning(check_seeds_in_catchments(s, catchments_sample), "Not all seeds are in catchments_sf")
})

test_that("check_seeds_in_catchments errors when no seeds are in catchments", {
  s <- tibble::tibble(CATCHNUM = c(1L, 2L), Areatarget = 1e9)
  expect_error(check_seeds_in_catchments(s, catchments_sample), "None of the seeds are in catchments_sf")
})


# fetch_builder_output
# write a fake BUILDER output table in a temporary folder
write_builder_csv <- function(dir, file, tab = builder_table_sample) {
  tab$OID <- seq_len(nrow(tab)) # BUILDER tables have an OID column
  utils::write.csv(tab, file.path(dir, file), row.names = FALSE)
}

test_that("fetch_builder_output reads the BENCHMARKS table and drops OID", {
  d <- withr::local_tempdir()
  write_builder_csv(d, "run1_COLUMN_All_Unique_BAs.csv")
  expect_message(fetch_builder_output(d), "run1_COLUMN_All_Unique_BAs.csv")
  out <- suppressMessages(fetch_builder_output(d))
  expect_s3_class(out, "tbl_df")
  expect_false("OID" %in% names(out))
  expect_equal(out, builder_table_sample, ignore_attr = TRUE)
})

test_that("fetch_builder_output reads UPSTREAM and DOWNSTREAM tables", {
  d <- withr::local_tempdir()
  write_builder_csv(d, "run1_UPSTREAM_CATCHMENTS_COLUMN.csv")
  write_builder_csv(d, "run1_DOWNSTREAM_CATCHMENTS_COLUMN.csv")
  expect_message(fetch_builder_output(d, "UPSTREAM"), "UPSTREAM_CATCHMENTS_COLUMN")
  expect_message(fetch_builder_output(d, "DOWNSTREAM"), "DOWNSTREAM_CATCHMENTS_COLUMN")
})

test_that("fetch_builder_output warns and returns the last table when several exist", {
  d <- withr::local_tempdir()
  write_builder_csv(d, "run1_COLUMN_All_Unique_BAs.csv")
  write_builder_csv(d, "run2_COLUMN_All_Unique_BAs.csv")
  expect_warning(
    expect_message(fetch_builder_output(d), "run2_COLUMN_All_Unique_BAs.csv"),
    "Multiple output tables"
  )
})

test_that("fetch_builder_output warns and returns NULL for an empty table", {
  d <- withr::local_tempdir()
  write_builder_csv(d, "run1_COLUMN_All_Unique_BAs.csv", builder_table_sample[0, ])
  expect_warning(out <- suppressMessages(fetch_builder_output(d)), "No conservation areas")
  expect_null(out)
})

test_that("fetch_builder_output ignores files that only contain the table name", {
  d <- withr::local_tempdir()
  write_builder_csv(d, "run1_COLUMN_All_Unique_BAs.csv")
  write_builder_csv(d, "run2_COLUMN_All_Unique_BAs.csv.bak")  # wrong extension
  write_builder_csv(d, "run3_COLUMN_All_Unique_BAsXcsv")     # "." must be literal
  expect_no_warning(expect_message(fetch_builder_output(d), "run1_COLUMN_All_Unique_BAs.csv"))
})

test_that("fetch_builder_output errors on bad inputs", {
  d <- withr::local_tempdir()
  expect_error(fetch_builder_output(NULL), "Builder_output directory")
  expect_error(fetch_builder_output(d, "OTHER"), "BENCHMARKS, UPSTREAM or DOWNSTREAM")
  expect_error(fetch_builder_output(d), "No output produced") # empty folder
})


# builder
# These tests stop before BUILDER is called, so they run without BenchmarkBuilder_cmd.exe
test_that("builder errors on invalid data_source", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  expect_error(
    builder(catchments_sample, data_source = "reserves", seeds = s, neighbours = NULL,
            builder_local_path = "not_used"),
    "should be one of"
  )
})

test_that("builder deletes its temp folder when it stops on an error", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  s$Areatarget[1] <- 0L # fails check_seeds_areatargets
  expect_error(
    suppressMessages(builder(catchments_sample, seeds = s, neighbours = NULL, builder_local_path = "not_used")),
    "All Areatarget values must be > 0"
  )
  expect_false(dir.exists(file.path(tempdir(), "builder")))
})

test_that("builder keeps out_dir when it stops on an error", {
  out <- file.path(withr::local_tempdir(), "builder_out")
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  s$Areatarget[1] <- 0L
  expect_error(
    suppressMessages(builder(catchments_sample, seeds = s, neighbours = NULL, out_dir = out, builder_local_path = "not_used")),
    "All Areatarget values must be > 0"
  )
  expect_true(dir.exists(out))
})

test_that("builder errors when out_dir has spaces", {
  s <- seeds(catchments_sample, areatarget_value = 1e9)
  expect_error(
    builder(catchments_sample, seeds = s, neighbours = NULL, out_dir = "C:/builder out", builder_local_path = "not_used"),
    "out_dir cannot have spaces"
  )
})
