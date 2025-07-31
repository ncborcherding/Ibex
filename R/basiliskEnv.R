#' @import basilisk
.get_ibex_env <- local({
  env <- NULL
  function() {
    if (is.null(env)) {
      env <<- basilisk::BasiliskEnvironment(
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
    }
    env
  }
})