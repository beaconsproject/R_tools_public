test_that("KS works", {
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
  expect_equal(calc_dissimilarity(reserves, "network", ref_poly, led_sample, 'continuous'),
               setNames(c(0.052, 0.137, 0.088), c("PB_0001", "PB_0002", "PB_0003")))
})

test_that("BC  works", {
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002", "PB_0003"))
  expect_equal(calc_dissimilarity(reserves, "network", ref_poly, led_sample, 'categorical', categorical_class_values = c(1,2,3), plot_out_dir = NULL, categorical_class_labels = c("low", "med", "high")),
               setNames(c(0.064279455, 0.079467331, 0.012132321), c("PB_0001", "PB_0002", "PB_0003")))
})


