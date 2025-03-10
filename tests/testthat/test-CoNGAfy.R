# test script for CoNGAfy.R - testcases are NOT comprehensive!

test_that("CoNGAfy works with Seurat object", {
  result <- CoNGAfy(ibex_example, method = "mean")
  
  expect_true(inherits(result, "SingleCellExperiment"))
  expect_gt(ncol(result), 0)
  expect_gt(nrow(result), 0)
})


test_that("CoNGAfy works with dist method", {
  result <- CoNGAfy(ibex_example, method = "dist")
  
  expect_true(inherits(result, "SingleCellExperiment"))
  expect_gt(ncol(result), 0)
  expect_gt(nrow(result), 0)
})

test_that("CoNGAfy filters cells correctly", {
  result <- CoNGAfy(ibex_example, method = "mean")
  expect_equal(ncol(result), 52) 
})

test_that("CoNGAfy stops if amino acid sequences are missing", {
  sc_example <- suppressWarnings(CreateSeuratObject(counts = matrix(rnorm(1000), 
                                                    nrow = 10, 
                                                    ncol = 100)))
  
  expect_error(CoNGAfy(sc_example, method = "mean"),
               "'CTaa' not found in this Seurat object\n ")
})

test_that("CoNGA.dist selects representative cells correctly", {
  result <- .CoNGA.dist(ibex_example, features = NULL, assay = "RNA")
  
  expect_true(inherits(result, "dgCMatrix"))
  expect_gt(ncol(result), 0)
  expect_gt(nrow(result), 0)
})

test_that("CoNGA.mean computes mean expression per clonotype", {
  result <- .CoNGA.mean(ibex_example, features = NULL, assay = "RNA")
  
  expect_true(inherits(result, "dgCMatrix"))
  expect_gt(ncol(result), 0)
  expect_gt(nrow(result), 0)
})
