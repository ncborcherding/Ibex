# test script for Ibex.matrix.R - testcases are NOT comprehensive!

test_that("maTrex works", {
  data("ibex_example")
  
  set.seed(42)
  
  ibex.result1 <- Ibex.matrix(ibex_example, 
                              chains = "Heavy",
                              method = "encoder",
                              encoder.model = "VAE",
                              encoder.input = "AF")
  
  expect_equal(
    ibex.result1,
    getdata("runIbex", "ibex.matrix_Heavy_VAE_AF")
  )
  
  ibex.result2 <- Ibex.matrix(ibex_example, 
                              chains = "Light",
                              method = "encoder",
                              encoder.model = "AE",
                              encoder.input = "OHE")
  
  expect_equal(
    ibex.result2,
    getdata("runIbex", "ibex.matrix_Light_AE_OHE")
  )
})