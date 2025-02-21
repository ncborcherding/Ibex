# test script for runTrex.R - testcases are NOT comprehensive!
library(Seurat)

test_that("runIbex handles incorrect inputs gracefully", {
  expect_error(runIbex(sc.data = ibex_example, chain = "Middle", method = "encoder"),
               "'arg' should be one of \"Heavy\", \"Light\"")
  expect_error(runIbex(sc.data = ibex_example, chain = "Heavy", method = "xyz"),
               "'arg' should be one of \"encoder\", \"geometric\"")
  expect_error(runIbex(sc.data = ibex_example, chain = "Heavy", method = "encoder", encoder.model = "ABC"),
               "'arg' should be one of \"CNN\", \"VAE\", \"CNN.EXP\", \"VAE.EXP\"")
  expect_error(runIbex(sc.data = ibex_example, chain = "Heavy", method = "encoder", encoder.input = "XYZ"),
               "arg' should be one of \"atchleyFactors\", \"crucianiProperties\", \"kideraFactors\", \"MSWHIM\", \"tScales\", \"OHE\"")
  expect_error(runIbex(sc.data = ibex_example, chain = "Heavy", method = "geometric", geometric.theta = "not_numeric"),
               "non-numeric argument to mathematical function")
})

keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping runIbex tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
    
  test_that("runIbex works with Seurat object", {
    suppressWarnings(sc_example <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10, ncol = 100)))
    sc_example[["CTaa"]] <- sample(c("CASSL", "CASST", NA, "NA_IGHV1", "None_IGHV2"), 100, replace = TRUE)
    sc_example[["CTgene"]] <- sample(c("NA_IGHV1.IGD1.IGJ1.IGM", "NA_IGHV1.IGD1.IGJ1.IGM", NA, "NA_IGHV1.IGD1.IGJ1.IGM", "None_IGHV1.IGD1.IGJ1.IGM"), 100, replace = TRUE)
    
    result <- runIbex(sc_example, 
                      chain = "Heavy", 
                      method = "encoder",
                      encoder.model = "VAE", 
                      encoder.input = "atchleyFactors",
                      reduction.name = "IbexTest", 
                      verbose = FALSE)
    
    expect_true("IbexTest" %in% names(result@reductions))
    expect_true(inherits(result, "Seurat"))
  })
  
  test_that("runIbex works with geometric method", {
    sc_example <- suppressWarnings(CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10, ncol = 100)))
    sc_example[["CTaa"]] <- sample(c("CASSL", "CASST", NA, "NA_IGHV1", "None_IGHV2"), 100, replace = TRUE)
    sc_example[["CTgene"]] <- sample(c("NA_IGHV1.IGD1.IGJ1.IGM", "NA_IGHV1.IGD1.IGJ1.IGM", NA, "NA_IGHV1.IGD1.IGJ1.IGM", "None_IGHV1.IGD1.IGJ1.IGM"), 100, replace = TRUE)
    
    result <- runIbex(sc_example, 
                      chain = "Heavy", 
                      method = "geometric",
                      geometric.theta = pi / 4, 
                      reduction.name = "IbexGeo", 
                      verbose = FALSE)
    
    expect_true("IbexGeo" %in% names(result@reductions))
    expect_true(inherits(result, "Seurat"))
  })
  
  test_that("runIbex filters cells correctly", {
    sc_example <- suppressWarnings(CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10, ncol = 100)))
    sc_example[["CTaa"]] <- c(rep("CASSL", 50), rep(NA, 50))
    sc_example[["CTgene"]] <- sample(c("NA_IGHV1.IGD1.IGJ1.IGM", "NA_IGHV1.IGD1.IGJ1.IGM", NA, "NA_IGHV1.IGD1.IGJ1.IGM", "None_IGHV1.IGD1.IGJ1.IGM"), 100, replace = TRUE)
    result <- runIbex(sc_example, 
                      chain = "Heavy", 
                      method = "encoder",
                      encoder.model = "VAE", 
                      encoder.input = "atchleyFactors",
                      reduction.name = "IbexFiltered", 
                      verbose = FALSE)
    
    expect_true("IbexFiltered" %in% names(result@reductions))
    expect_lt(ncol(result), 100)  # Ensures some cells were filtered out
  })
  
  test_that("runIbex stops if amino acid sequences are missing", {
    sc_example <- suppressWarnings(CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10, ncol = 100)))
    
    expect_error(runIbex(sc_example, 
                         chain = "Heavy", 
                         method = "encoder",
                         encoder.model = "VAE", 
                         encoder.input = "atchleyFactors", 
                         verbose = FALSE),
                 "Amino acid sequences are not added to the single-cell object correctly.")
  })
  
  test_that("runIbex works with different reduction names", {
    sc_example <- suppressWarnings(CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10, ncol = 100)))
    sc_example[["CTaa"]] <- sample(c("CASSL", "CASST", NA, "NA_IGHV1", "None_IGHV2"), 100, replace = TRUE)
    sc_example[["CTgene"]] <- sample(c("NA_IGHV1.IGD1.IGJ1.IGM", "NA_IGHV1.IGD1.IGJ1.IGM", NA, "NA_IGHV1.IGD1.IGJ1.IGM", "None_IGHV1.IGD1.IGJ1.IGM"), 100, replace = TRUE)
    result1 <- runIbex(sc_example, 
                       chain = "Heavy", 
                       method = "encoder",
                       encoder.model = "VAE", 
                       encoder.input = "atchleyFactors",
                       reduction.name = "Ibex1", 
                       verbose = FALSE)
    
    result2 <- runIbex(sc_example, chain = "Heavy", 
                       method = "encoder",
                       encoder.model = "VAE", 
                       encoder.input = "atchleyFactors",
                       reduction.name = "Ibex2", 
                       verbose = FALSE)
    
    expect_true("Ibex1" %in% names(result1@reductions))
    expect_true("Ibex2" %in% names(result2@reductions))
  })
}
