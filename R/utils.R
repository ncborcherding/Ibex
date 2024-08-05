"%!in%" <- Negate("%in%")

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

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

#Shhhhhh
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#This is to check that all the cdr3 sequences are < 70 residues
checkLength <- function(x) {
  if(any(na.omit(nchar(x[,"cdr3_aa"])) > 70)) {
    stop("Models have been trained on cdr3 sequences 
         less than 70 amino acid residues. Please
         filter the larger sequences before running")
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

one.hot.organizer <- function(refer) {
  reference <- ibex.data[[1]]
  int <- matrix(ncol = length(reference$aa) + 1, nrow = length(refer))
  for(i in seq_along(refer)) {
    if (is.na(refer[i])) {
      next()
    }
    int[i,which(reference$aa %in% refer[i])] <- 1
  }
  int[is.na(int)] <- 0
  return(int)
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

AF.col <- c(2,3,4,5,6)
KF.col <- c(7,8,9,10,11,12,13,14,15,16)

#Generates the 30 vector based on autoencoder model 
auto.embedder <- function(array.reshape, aa.model, encoder.input) {
  array.reshape[is.na(array.reshape)] <- 0
  score <- stats::predict(aa.model, t(array.reshape), verbose = 0)
  return(score)
}

dist.convert <- function(dist_obj, k) {
  if (length(k) > 1) stop("The function is not 'vectorized'!")
  n <- attr(dist_obj, "Size")
  if (k < 1 || k > n) stop("k out of bound!")
  ##
  i <- 1:(k - 1)
  j <- rep.int(k, k - 1)
  v1 <- dist_obj[f(j, i, dist_obj)]
  return(v1)
}
