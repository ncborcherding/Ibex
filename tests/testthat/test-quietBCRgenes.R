# test script for quietBCRgenes.R - testcases are NOT comprehensive!

test_that("quietBCRgenes works", {
  
  data("ibex_example")
  
  features <- rownames(ibex_example@assays$RNA$counts)
  
  expect_equal(
    quietBCRgenes(features),
    getdata("quietBCRgenes", "quietBCRgenes_feature.vector")
  )
  
  ibex_example@assays$RNA@var.features <- features
  Seurat::DefaultAssay(ibex_example) <- "RNA"
  
  ibex_example <- quietBCRgenes(ibex_example)
  
  
  expect_equal(
    quietBCRgenes(features),
    ibex_example@assays$RNA@var.features
  )
})