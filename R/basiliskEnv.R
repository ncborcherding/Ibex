#' @import basilisk
IbexEnv <- BasiliskEnvironment(
  envname = "IbexEnv",
  pkgname = "Ibex",
  packages = c(
    "python=3.9",
    "tensorflow=2.13.*",   
    "h5py=3.13",
    "numpy=1.26")
)