.onLoad <- function(libname, pkgname) {
  # Set a safe, temporary cache directory for basilisk if one is not already set.
  # This is critical for automated build environments (like Bioconductor's)
  # where default user-level cache directories are not writable.
  if (nchar(Sys.getenv("BASILISK_CACHE_DIR")) == 0) {
    cache_dir <- file.path(tempdir(), "basilisk_cache")
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    Sys.setenv(BASILISK_CACHE_DIR = cache_dir)
  }
}