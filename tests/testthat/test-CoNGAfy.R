# test script for CoNGAfy.R - testcases are NOT comprehensive!

test_that("CoNGAfy works", {
  
  data("ibex_example")
  
  conga_reduction <- CoNGAfy(ibex_example)
  
  expect_equal(
    conga_reduction@meta.data,
    getdata("CoNGAfy", "CoNGAfy_meta.data")
  )
  
  expect_equal(
    conga_reduction@assays$RNA@layers$counts,
    getdata("CoNGAfy", "CoNGAfy_counts"),
    tolerance=1e-2
  )
  
  conga_mean_reduction <- CoNGAfy(ibex_example,
                                  method = "mean")
  
  expect_equal(
    conga_mean_reduction@meta.data,
    getdata("CoNGAfy", "CoNGAfy_mean_meta.data")
  )
  
  expect_equal(
    conga_mean_reduction@assays$RNA@layers$counts,
    getdata("CoNGAfy", "CoNGAfy_mean_counts"),
    tolerance=1e-2
  )
})