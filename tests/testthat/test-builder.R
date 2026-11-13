
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
