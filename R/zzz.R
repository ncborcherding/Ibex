.onLoad <- function(libname, pkgname) {
  # Create a stable cache dir for basilisk
  base <- tools::R_user_dir("Ibex", which = "cache")
  exdir <- file.path(base, "basilisk")
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  
  # Tell basilisk and dir.expiry to use it
  options(basilisk.external.dir = exdir)
  options(dir.expiry.dir = exdir)
}