library(sf)

test_that("projection error catches", {
  expect_error(
    criteria_to_catchments(st_transform(catchments_sample, 3857), led_sample, "led"),
    "is not TRUE"
    )
})

test_that("CATCHNUM error catches", {
  expect_error(
    criteria_to_catchments(catchments_sample[c(2:14)], led_sample, "led"),
    "CATCHNUM"
    )
})

test_that("areas are as expected", {
  expect_snapshot_output(as.data.frame(criteria_to_catchments(catchments_sample["CATCHNUM"], led_sample, "led", c(1:15))))
})

test_that("NAs are removed", {
  led_nas <- led_sample
  led_nas[led_nas == 5] <- NA
  expect_snapshot_output(as.data.frame(criteria_to_catchments(catchments_sample["CATCHNUM"], led_nas, "led", c(1:15))))
})

test_that("zeros are returned when NA is only value in catchment", {
  led_nas <- led_sample
  led_nas[led_nas > 0] <- NA
  expect_snapshot_output(as.data.frame(criteria_to_catchments(catchments_sample["CATCHNUM"], led_nas, "led", c(1:15))))
})

