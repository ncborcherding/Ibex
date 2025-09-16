# test script for Ibex_matrix.R - testcases are NOT comprehensive!
library(Ibex)
ibex_example <- get(data("ibex_example"))

test_that("Ibex_matrix handles incorrect inputs gracefully", {
  expect_error(Ibex_matrix(input.data = ibex_example, chain = "Middle", method = "encoder"),
               "'arg' should be one of \"Heavy\", \"Light\"")
  expect_error(Ibex_matrix(input.data = ibex_example, chain = "Heavy", method = "xyz"),
               "'arg' should be one of \"encoder\", \"geometric\"")
  expect_error(Ibex_matrix(input.data = ibex_example, chain = "Heavy", method = "encoder", encoder.model = "ABC"),
               "'arg' should be one of \"CNN\", \"VAE\", \"CNN.EXP\", \"VAE.EXP\"")
  expect_error(Ibex_matrix(input.data = ibex_example, chain = "Heavy", method = "encoder", encoder.input = "XYZ"),
               "arg' should be one of \"atchleyFactors\", \"crucianiProperties\", \"kideraFactors\", \"MSWHIM\", \"tScales\", \"OHE\"")
  expect_error(Ibex_matrix(input.data = ibex_example, chain = "Heavy", method = "geometric", geometric.theta = "not_numeric"),
               "non-numeric argument to mathematical function")
})

keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping Ibex_matrix tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
    
  test_that("Ibex_matrix returns expected output format", {
    result <- Ibex_matrix(input.data = ibex_example, 
                          chain = "Heavy", 
                          method = "encoder",
                          encoder.model = "VAE", 
                          encoder.input = "atchleyFactors", 
                          verbose = FALSE)
    expect_true(is.data.frame(result))
    expect_true(all(grepl("^Ibex_", colnames(result))))
    expect_gt(nrow(result), 0)
    expect_gt(ncol(result), 0)
  })
  
  test_that("Ibex_matrix works with encoder method", {
    result <- Ibex_matrix(input.data = ibex_example, 
                          chain = "Light", 
                          method = "encoder",
                          encoder.model = "CNN", 
                          encoder.input = "OHE", 
                          verbose = FALSE)
    expect_true(is.data.frame(result))
    expect_true(all(grepl("^Ibex_", colnames(result))))
  })
  
  test_that("Ibex_matrix works with geometric method", {
    result <- Ibex_matrix(input.data = ibex_example, 
                          chain = "Heavy", 
                          method = "geometric",
                          geometric.theta = pi / 4, 
                          verbose = FALSE)
    expect_true(is.data.frame(result))
    expect_true(all(grepl("^Ibex_", colnames(result))))
  })
  
  test_that("Ibex_matrix handles different species options", {
    result1 <- Ibex_matrix(input.data = ibex_example, 
                           chain = "Heavy", 
                           method = "encoder",
                           encoder.model = "VAE", 
                           encoder.input = "atchleyFactors", 
                           species = "Human", 
                           verbose = FALSE)
    result2 <- Ibex_matrix(input.data = ibex_example, 
                           chain = "Heavy", 
                           method = "encoder",
                           encoder.model = "VAE", 
                           encoder.input = "atchleyFactors", 
                           species = "Mouse", 
                           verbose = FALSE)
    expect_true(is.data.frame(result1))
    expect_true(is.data.frame(result2))
    expect_true(all(grepl("^Ibex_", colnames(result1))))
    expect_true(all(grepl("^Ibex_", colnames(result2))))
  })
}
