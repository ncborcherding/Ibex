# test script for runTrex.R - testcases are NOT comprehensive!

test_that("runIbex works with seurat objects", {
	data("ibex_example")
  set.seed(42)
  
  ibex_example <- runIbex(ibex_example, 
                          chains = "Heavy",
                          method = "encoder",
                          encoder.model = "VAE",
                          encoder.input = "AF", 
                          reduction.name = "Heavy_VAE_AF")
                       
	expect_equal(
		ibex_example@reductions$Heavy_VAE_AF@cell.embeddings,
		getdata("runIbex", "runIbex_Heavy_VAE_AF_reduction")
	)
	
	ibex_example <- runIbex(ibex_example, 
	                        chains = "Light",
	                        method = "encoder",
	                        encoder.model = "AE",
	                        encoder.input = "KF", 
	                        reduction.name = "Light_AE_KF")
	
	expect_equal(
	  ibex_example@reductions$Light_AE_KF@cell.embeddings,
	  getdata("runIbex", "runIbex_Light_AE_KF_reduction")
	)
	
	ibex_example <- runIbex(ibex_example, 
	                        chains = "Heavy",
	                        method = "encoder",
	                        encoder.model = "VAE",
	                        encoder.input = "OHE", 
	                        reduction.name = "Heavy_VAE_OHE")
	
	expect_equal(
	  ibex_example@reductions$Heavy_VAE_OHE@cell.embeddings,
	  getdata("runIbex", "runIbex_Heavy_VAE_OHE_reduction")
	)
	
})
