# test script for combineExpandedBCR.R - testcases are NOT comprehensive!
library(Ibex)
ibex_vdj <- get(data("ibex_vdj"))

test_that("combineExpandedBCR handles incorrect input gracefully", {
  expect_error(combineExpandedBCR(NULL, samples = "Sample1"),
               "Input data must be a list of data frames.")
  
  invalid_data <- list(data.frame(cdr1 = c("AA", "BB"), cdr3 = c("CC", "DD")))
  expect_error(combineExpandedBCR(invalid_data, samples = "Sample1"),
               "Each data frame must contain 'cdr1', 'cdr2', and 'cdr3' columns.")
})

test_that("combineExpandedBCR correctly concatenates CDR sequences", {
  
  modified_data <- combineExpandedBCR(list(ibex_vdj), samples = "Sample1")
  
  expect_true(any(grepl("-", modified_data[[1]]$CTaa)))
})

test_that("combineExpandedBCR integrates correctly with combineBCR", {
  
  result <- combineExpandedBCR(list(ibex_vdj), samples = "Sample1")
  expect_true(is.list(result))
  expect_true(all(c("barcode", "CTaa") %in% colnames(result[[1]])))
  expect_gt(nrow(result[[1]]), 0)
})

test_that("combineExpandedBCR correctly assigns sample labels", {
  
  result <- combineExpandedBCR(list(ibex_vdj), samples = "Sample1")
  
  expect_true("sample" %in% colnames(result[[1]]))
  expect_equal(result[[1]]$sample[1], "Sample1")
})

test_that("combineExpandedBCR handles multiple sample inputs correctly", {
  
  result <- combineExpandedBCR(list(ibex_vdj, ibex_vdj), samples = c("Sample1", "Sample2"))
  
  expect_true(length(result) == 2)
  expect_equal(result[[1]]$sample[1], "Sample1")
  expect_equal(result[[2]]$sample[1], "Sample2")
})


