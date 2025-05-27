#' @import basilisk
IbexEnv <- BasiliskEnvironment(
  envname = "IbexEnv",
  pkgname = "Ibex",
  packages = c(
    "python=3.9",
    "keras=3.6.*",            
    "tensorflow=2.18.*",      
    "h5py=3.13",
    "numpy=1.26"
  )
)