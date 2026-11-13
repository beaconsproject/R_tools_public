
# gen_network_names
test_that("k of 2 works", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3"), 2),
    c("PB_1__PB_2", "PB_1__PB_3", "PB_2__PB_3")
    )
})

test_that("k of 3 works", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_4"), 3),
    c("PB_1__PB_2__PB_3", "PB_1__PB_2__PB_4", "PB_1__PB_3__PB_4", "PB_2__PB_3__PB_4")
    )
})

test_that("sort is working", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_4"), 3),
    gen_network_names(c("PB_4", "PB_3", "PB_2", "PB_1"), 3)
    )
})

test_that("force_in is working with one value", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2, "PB_1"),
    c("PB_1__PB_11", "PB_1__PB_2", "PB_1__PB_3")
    )
})

test_that("force_in is working with vector", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2, c("PB_1", "PB_2")),
    c("PB_1__PB_11", "PB_1__PB_2", "PB_1__PB_3", "PB_2__PB_3", "PB_11__PB_2")
    )
})

test_that("force_in is working with multi-name string", {
  expect_equal(
    gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 3, c("PB_1__PB_11")),
    c("PB_1__PB_11__PB_2", "PB_1__PB_11__PB_3")
    )
})


# sep_network_names
test_that("result for string", {
  expect_equal(
    sep_network_names("PB_1__PB_2"),
    c("PB_1", "PB_2")
  )
})

test_that("result for vector", {
  expect_equal(
    sep_network_names(c("PB_1__PB_2", "PB_1__PB_3")),
    list("PB_1__PB_2" = c("PB_1", "PB_2"), "PB_1__PB_3" = c("PB_1", "PB_3"))
  )
})