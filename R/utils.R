"%!in%" <- Negate("%in%")

amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

# Add to meta data some of the metrics calculated
#' @importFrom rlang %||%
#' @importFrom SingleCellExperiment colData
add.meta.data <- function(sc, meta, header) {
if (inherits(x=sc, what ="Seurat")) { 
  col.name <- names(meta) %||% colnames(meta)
  sc[[col.name]] <- meta
} else {
  rownames <- rownames(colData(sc))
  colData(sc) <- cbind(colData(sc), 
          meta[rownames,])[, union(colnames(colData(sc)),  colnames(meta))]
  rownames(colData(sc)) <- rownames  
}
  return(sc)
}

# This is to grab the metadata from a Seurat or SCE object
#' @importFrom SingleCellExperiment colData 
#' @importFrom methods slot
grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="SingleCellExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}

# This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

# This is to check that all the CDR3 sequences are < 45 residues or < 90 for CDR1/2/3
#' @importFrom stats na.omit
checkLength <- function(x, expanded = NULL) {
  cutoff <- ifelse( expanded == FALSE || is.null(expanded), 45, 90)
  if(any(na.omit(nchar(x)) > cutoff)) {
    stop("Models have been trained on sequences less than ", cutoff, 
         " amino acid residues. Please filter the larger sequences before running")
  }
}
# Returns appropriate encoder model
#' @importFrom utils download.file read.csv
#' @importFrom tools R_user_dir
#' @importFrom utils download.file read.csv
#' @importFrom tools R_user_dir
aa.model.loader <- function(species,
                            chain,
                            encoder.input,
                            encoder.model) {

  ## 1. Expected filename
  file_name <- paste0(
    species, "_", chain, "_",
    encoder.model, "_", encoder.input,
    "_encoder.keras")
  
  ## 2. Sanity-check against metadata.csv 
  meta <- read.csv(
    system.file("extdata", "metadata.csv", package = "Ibex"),
    stringsAsFactors = FALSE
  )
  
  if (!file_name %in% meta[[1]])
    stop("Model '", file_name, "' is not listed in metadata.csv.")
  
  ## 3. Cache directory 
  cache_dir <- tools::R_user_dir("Ibex", which = "cache")
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  local_path <- file.path(cache_dir, file_name)
  
  ## 4. Download if we have never seen this model before
  if (!file.exists(local_path)) {
    message("Downloading model '", file_name, " ...")
    base_url <- "https://zenodo.org/record/14919286/files"
    status   <- utils::download.file(
      url      = file.path(base_url, file_name),
      destfile = local_path,
      mode     = "wb",
      quiet    = TRUE
    )
    if (status != 0)
      stop("Download of model '", file_name,
           "' failed (status ", status, ").")
  }
  
  ## 5. Done return the path for use in basiliskRun()
  normalizePath(local_path, winslash = "/")
}



# Add the dimRed to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
adding.DR <- function(sc, reduction, reduction.name) {
  if (inherits(sc, "Seurat")) {
    DR <- suppressWarnings(CreateDimReducObject(
      embeddings = as.matrix(reduction),
      loadings = as.matrix(reduction),
      projected = as.matrix(reduction),
      stdev = rep(0, ncol(reduction)),
      key = reduction.name,
      jackstraw = NULL,
      misc = list()))
    sc[[reduction.name]] <- DR
  } else if (inherits(sc, "SingleCellExperiment")) {
    reducedDim(sc, reduction.name) <- reduction
  }
  return(sc)
}

# Get the external dir in a way that won't explode during lazy-load on older builders.
.ibex_get_external_dir <- function() {
  if (requireNamespace("basilisk.utils", quietly = TRUE)) {
    ns <- asNamespace("basilisk.utils")
    if (exists("getExternalDir", envir = ns, inherits = FALSE)) {
      return(get("getExternalDir", envir = ns)())
    }
  }
  # Fallback that works on the build farm without env vars.
  file.path(tempdir(), "Ibex", "basilisk")
}

# Ensure the directory exists & is writable; return logical.
.ibex_ensure_external_dir <- function() {
  exdir <- .ibex_get_external_dir()
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  file.access(exdir, 2L) == 0L
}
