# avg_rank
test_that("ranks are as expected", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(10, 1, 8, 11, 1))
  
  expect_equal(
    avg_rank(tbl, 'col1', 'col2'),
    c(2, 2, 3.5, 4.5, 3)
  )
})

test_that("test all NA error", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(NA, NA, NA, NA, NA))
  
  expect_error(
    avg_rank(tbl, c('col1', 'col3'), 'col2'),
    "All values in column"
  )
})

test_that("test some NA", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(NA, NA, 1, NA, NA))
  
  expect_warning(
    expect_equal(
      avg_rank(tbl, c('col1', 'col3'), 'col2'),
      c(2.5, 2.5, 2.67, 4.17, 3.17)
    ),
    "contains NA values"
  )
})

test_that("test single column asc", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(NA, NA, 1, NA, NA))
  
  expect_equal(
    avg_rank(tbl, columns_asc = 'col1'),
    c(1,2,3,4,5)
  )
})

test_that("test single column desc", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(NA, NA, 1, NA, NA))
  
  expect_equal(
    avg_rank(tbl, columns_desc = 'col1'),
    c(5,4,3,2,1)
  )
})

test_that("test wrong column name", {
  tbl <- data.frame(col1 = c(1,2,3,4,5), col2 = c(0.4, 0.5, 0.2, 0.05, 0.9), col3 = c(NA, NA, 1, NA, NA))
  
  expect_error(
    avg_rank(tbl, columns_desc = 'col4'),
    "Columns not in table"
  )
})
