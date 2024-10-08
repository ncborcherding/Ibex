#' @param input.data Single Cell Object in Seurat or SingleCell Experiment format or
#' the output of combineBCR() in scRepertoire
#' @param chains The immune receptor chain to encode: "Heavy", "Light", or "Both"
#' @param encoder.function The method to prepare the sequencing information - 
#' "onehotEncoder", "propertyEncoder", or "geometricEncoder"
#' @param aa.method.to.use The method or approach to use for the conversion:
#' \itemize{
#'   \item{Individual sets: atchleyFactors, crucianiProperties, FASGAI, kideraFactors, MSWHIM,
#'   ProtFP, stScales, tScales, VHSE, zScales"}
#'   \item{Multiple Sets: c("atchleyFactors", "VHSE") }
#' } 
#' 
#' @importFrom immApex getIR onehotEncoder propertyEncoder geometricEncoder

Ibex.matrix <- function(input.data, 
                        chains = "Heavy", 
                        encoder.function = "onehotEncoder",
                        aa.method.to.use = NULL,
                        encoder.model = "VAE",
                        geometric.matrix = "BLOSUM62",
                        theta = pi/3) {
  
  chains <- chain.checker(chains)
  if(chains %in% c("Heavy", "Light")) {
    BCR <- getIR(input.data, chains)
  } else if(chains == "Both") {
    lapply(c("Heavy", "Light"), function(x) {
      tmp <- getIR(input.data, x)[[1]]
    }) -> BCR
    names(BCR) <- c("Heavy", "Light")
  }
  #TODO Extract Sequences/Combine Sequences
  checkLength(BCR[[1]])
  
  if(encoder.function %in% c("onehotEncoder", "propertyEncoder")) {                       
    sequence.matrix <- switch(encoder.function,
                              "onehotEncoder" = onehotEncoder(input.sequences, 
                                                              sequence.dictionary = sequence.dictionary,
                                                              convert.to.matrix = TRUE),
                              "propertyEncoder" = propertyEncoder(input.sequences, 
                                                                  method.to.use = aa.method.to.use,
                                                                  convert.to.matrix = TRUE))
    
    aa.model <- quiet(aa.model.loader(chains, encoder.function, encoder.model))
    membership <- BCR[[1]]
    names <- membership$barcode
    sequence.matrix <- auto.embedder(sequence.matrix, aa.model, encoder.input)
    
    sequence.matrix <- data.frame(sequence.matrix)
    rownames(sequence.matrix) <- unique(membership[,"barcode"])
                      
  } else if(encoder.function == "geometricEncoder") {
    sequence.matrix <- geometricEncoder(input.sequences, 
                                        theta  = theta,
                                        method.to.use = geometric.matrix,
                                        convert.to.matrix = TRUE)
  }
  colnames(sequence.matrix) <- paste0("Ibex_", seq_len(ncol(sequence.matrix)))
  return(sequence.matrix)
}


#' @importFrom SingleCellExperiment colData 
#' @importFrom methods slot
.grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
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


