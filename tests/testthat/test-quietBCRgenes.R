# test script for quietBCRgenes.R - testcases are NOT comprehensive!

test_that("quietBCRgenes works", {
  
  data("ibex_example")
  
  features <- rownames(ibex_example@assays$RNA$counts)
  
  expect_equal(
    quietBCRgenes(features),
    getdata("quietBCRgenes", "quietBCRgenes_feature.vector")
  )
  
  SeuratObject::DefaultAssay(ibex_example) <- "RNA"
  SeuratObject::VariableFeatures(ibex_example, assay="RNA") <- features
  
  ibex_example <- quietBCRgenes(ibex_example)
  
  expect_equal(
    quietBCRgenes(features),
    Seurat::VariableFeatures(ibex_example)
  )
})