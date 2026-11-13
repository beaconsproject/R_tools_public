test_that("benchmark grouping works", {
  reserves <- dissolve_catchments_from_table(
    catchments_sample, 
    builder_table_sample,
    "network")
  
  expect_equal(group_conservation_areas(reserves, 10000),
               c(1,2,3)
               )
})

test_that("benchmark grouping works with bigger grid", {
  reserves <- dissolve_catchments_from_table(
    catchments_sample, 
    builder_table_sample,
    "network")
  
  expect_equal(group_conservation_areas(reserves, 100000),
               c(1,1,1)
               )
})


test_that("network grouping works", {
  reserves <- dissolve_catchments_from_table(
   catchments_sample, 
   builder_table_sample,
   "network")
  network_names <- gen_network_names(reserves$network, 2)
  networks_sf <- build_network_polygons(reserves, "network", network_names)
  
  expect_equal(group_conservation_area_networks(networks_sf$network, reserves, "network", 10000),
               c(1,2,3)
               )
})

test_that("network grouping works with bigger grid", {
  reserves <- dissolve_catchments_from_table(
   catchments_sample, 
   builder_table_sample,
   "network")
  network_names <- gen_network_names(reserves$network, 2)
  networks_sf <- build_network_polygons(reserves, "network", network_names)
  
  expect_equal(group_conservation_area_networks(networks_sf$network, reserves, "network", 100000),
               c(1,1,1)
               )
})