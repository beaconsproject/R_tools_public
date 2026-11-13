
# gen_targets
test_that("targets are as expected with 5000 reserve size", {
  expect_snapshot(gen_targets(ref_poly, led_sample, 5000))
})

test_that("targets are as expected with 10000 reserve size", {
  expect_snapshot(gen_targets(ref_poly, led_sample, 10000))
})

test_that("targets are as expected with class_vals", {
  expect_snapshot(gen_targets(ref_poly, led_sample, 10000, c(1:5)))
})



# evaluate_targets_using_catchments
test_that("evaluation table is as expected", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_snapshot(as.data.frame(evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, colnames(builder_table_sample))))
})

test_that("unique catchments work in evaluation", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_equal(
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001"))$prop_target_met,
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001__PB_0001"))$prop_target_met
    )
})

test_that("benchmark order does not matter", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_equal(
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001__PB_0002"))$prop_target_met,
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0002__PB_0001"))$prop_target_met
    )
})

test_that("default network_names works", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_equal(
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001","PB_0002","PB_0003")),
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table)
    )
})

test_that("warning is given when all target classes are not in catchments", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  catchments_test <- catchments_sample[names(catchments_sample)[names(catchments_sample) != "led_10"]]
  expect_warning(
    evaluate_targets_using_catchments(catchments_test, "led", builder_table_sample, target_table),
    "Classes are not in catchments: led_10"
    )
})

test_that("OID is removed", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  builder_table_sample_oid <- builder_table_sample %>%
    dplyr::mutate(OID = 1:nrow(builder_table_sample))
  
  expect_equal(
    as.data.frame(evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table)),
    as.data.frame(evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample_oid, target_table))
    )
})

# evaluate_targets_using_clip
test_that("evaluation table is as expected", {
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002"))
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_snapshot(as.data.frame(evaluate_targets_using_clip(reserves, "network", led_sample, target_table)))
})

test_that("benchmarks result matches evaluate_targets_using_catchments()", {
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002"))
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_equal(
    evaluate_targets_using_clip(reserves, "network", led_sample, target_table),
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001","PB_0002")),
    tolerance = 1e-5
    )
})

test_that("networks result matches evaluate_targets_using_catchments()", {
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001__PB_0002"))
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  expect_equal(
    evaluate_targets_using_clip(reserves, "network", led_sample, target_table),
    evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, c("PB_0001__PB_0002")),
    tolerance = 1e-5
    )
})

# evaluate_targets_using_networks
test_that("evaluation tables are as expected for 2 benchmark networks", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  benchmarks_results <- evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, colnames(builder_table_sample))
  network_list <- gen_network_names(unique(benchmarks_results$network), 2)
  expect_snapshot(as.data.frame(evaluate_targets_using_networks(benchmarks_results, network_list)))
})

test_that("evaluation tables are as expected for 3 benchmark network", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  benchmarks_results <- evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, colnames(builder_table_sample))
  expect_snapshot(evaluate_targets_using_networks(benchmarks_results, "PB_0001__PB_0002__PB_0003"))
})

test_that("network error catches", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  benchmarks_results <- evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, colnames(builder_table_sample))
  names(benchmarks_results)[names(benchmarks_results)=="network"] <- "networks"
  expect_error(
    evaluate_targets_using_networks(benchmarks_results, "PB_0001__PB_0002__PB_0003"),
    "network"
  )
})

test_that("eval table error catches", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  benchmarks_results <- evaluate_targets_using_catchments(catchments_sample, "led", builder_table_sample, target_table, colnames(builder_table_sample))
  names(benchmarks_results)[names(benchmarks_results)=="class_value"] <- "class"
  expect_error(
    evaluate_targets_using_networks(benchmarks_results, "PB_0001__PB_0002__PB_0003"),
    "network_evaluation_table must contain"
  )
})


# summarize_representation_results
test_that("summary gaps are correct", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002"))
  network_evaluation_table <- evaluate_targets_using_clip(reserves, "network", led_sample, target_table)
  expect_equal(
    list(
      summarize_representation_results(network_evaluation_table, "led")$led_gaps,
      summarize_representation_results(network_evaluation_table, "led", 0.9, 0.05)$led_gaps,
      summarize_representation_results(network_evaluation_table, "led", gaps = FALSE)$led_passed,
      summarize_representation_results(network_evaluation_table, "led", 0.9, 0.05, gaps = FALSE)$led_passed),
    list(
      c(12, 6),
      c(2, 1),
      c(3, 9),
      c(1, 2))
    )
})

test_that("summary names are correct", {
  target_table <- gen_targets(ref_poly, led_sample, 1600)
  reserves <- dissolve_catchments_from_table(catchments_sample, builder_table_sample, "network", dissolve_list = c("PB_0001", "PB_0002"))
  network_evaluation_table <- evaluate_targets_using_clip(reserves, "network", led_sample, target_table)
  expect_equal(
    c(colnames(summarize_representation_results(network_evaluation_table, "led"))[2],
      colnames(summarize_representation_results(network_evaluation_table, "led", 0.9, 0.05, suffix = "_test123"))[2],
      colnames(summarize_representation_results(network_evaluation_table, "led", gaps = FALSE))[2],
      colnames(summarize_representation_results(network_evaluation_table, "led", 0.9, 0.05, gaps = FALSE, suffix = "_test123"))[2]),
    c("led_gaps", "led_gaps_test123", "led_passed", "led_passed_test123")
  )
})
