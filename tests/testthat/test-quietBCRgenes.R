# test script for quietBCRgenes.R - testcases are NOT comprehensive!

test_that("quietBCRgenes works", {
  
  data("ibex_example")
  
  features <- rownames(ibex_example@assays@data$counts)
  
  expect_equal(
    quietBCRgenes(features),
    getdata("quietBCRgenes", "quietBCRgenes_feature.vector")
  )
})