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

#All the receptor chains shall be capitlized
#' @importFrom stringr str_to_title
chain.checker <- function(chain) {
  chain <- str_to_title(chain)
  if(chain %!in% c("Heavy", "Light", "Both")) {
    stop("Please select 'Heavy', 'Light' or 'Both' for the chains parameter,")
  }
}

#Function to pull and organize BCR depending on the chain selected
#' @importFrom stringr str_split
getBCR <- function(sc, chains) {
  if (inherits(x=sc, what ="Seurat") | inherits(x=sc, what ="SingleCellExperiment")) {
    meta <- grabMeta(sc)
  } else {
    meta <- do.call(rbind,sc)
    rownames(meta) <- meta[,"barcode"]
  }
  tmp <- data.frame(barcode = rownames(meta), 
                    str_split(meta[,"CTaa"], "_", simplify = TRUE), 
                    str_split(meta[,"CTgene"], "_", simplify = TRUE))
  if (length(chains) == 1 && chains != "both") {
    if (chains %in% c("Heavy")) { #here
      pos <- list(c(2,4))
    } else if (chains %in% c("Light")) { #here
      pos <- list(c(3,5))
    }
  } else {
    pos <- list(one = c(2,4), two = c(3,5))
    ch.1 <- grep("IGH|IGL",sc[[]]$CTgene[1]) #here
    chains <- c("heavy", "light") #here
  }
  BCR <- NULL
  for (i in seq_along(pos)) {
    sub <- as.data.frame(tmp[,c(1,pos[[i]])])
    colnames(sub) <- c("barcode", "cdr3_aa", "genes")
    sub$v <- str_split(sub$genes, "[.]", simplify = TRUE)[,1]
    sub$j <- str_split(sub$genes, "[.]", simplify = TRUE)[,2]
    sub[sub == ""] <- NA
    BCR[[i]] <- sub
    sub <- NULL
  }
  names(BCR) <- chains
  return(BCR)
}

#This is to grab the metadata from a Seurat or SCE object
#' @importFrom SingleCellExperiment colData 
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

#This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#This is to check that all the cdr3 sequences are < 45 residues or < 90 for cdr1/2/3
checkLength <- function(x, expanded = NULL) {
  cutoff <- ifelse(is.null(expanded)) | expanded == FALSE, 45, 90)
#TODO will need to update column pulling for expanded sequences
  if(any(na.omit(nchar(x[,"cdr3_aa"])) > cut.off)) {
    stop(paste0("Models have been trained on cdr3 sequences 
         less than ", cutoff, " amino acid residues. Please
         filter the larger sequences before running"))
  }
}
#Returns appropriate model for autoencoder
#' @importFrom tensorflow tf
#' @importFrom keras load_model_hdf5
aa.model.loader <- function(chain, encoder.input, encoder.model) {
    #TODO allow for multiple species when I get models trained
    species <- "Human"
    select  <- system.file("extdata", paste0(species, "_", chain, "_", 
                               encoder.input, "_", encoder.model, ".h5"), 
                          package = "Ibex")
    model <- quiet(load_model_hdf5(select, compile = FALSE))
    return(model)
}


#Add the eigenvectors to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
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

