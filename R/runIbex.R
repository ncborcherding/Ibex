#' Main Ibex interface
#' 
#' Use this to run the Ibex algorithm to return latent vectors in 
#' the form of a matrix or if you prefer.
#' 
#' @examples
#'ibex_values <- Ibex.matrix(ibex_example, 
#'                           chains = "Heavy",
#'                           method = "encoder",
#'                           encoder.model = "VAE",
#'                           encoder.input = "AF")
#'                           
#'ibex_values <- Ibex.matrix(ibex_example, 
#'                           chains = "Heavy",
#'                           method = "geometric",
#'                           theta = pi)
#'                         
#' @param input.data Single Cell Object in Seurat or Single Cell Experiment format or
#' the output of combineBCR() in scRepertoire
#' @param chain Heavy or Light
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on the BLOSUM62 matrix
#' @param encoder.model "CNN" = convolutional neural network-based autoencoder, 
#' "VAE" = variation autoencoder, or "biLSTM" for Bidirectional LSTM
#' @param encoder.input Type of model to use
#'   \itemize{
#'     \item{Amino Acid Properties: "atchleyFactors", "crucianiProperties", "FASGAI", "kideraFactors", "MSWHIM", "ProtFP", "stScales", "tScales", "VHSE", "zScales"}
#'     \item{"OHE" for One Hot Encoding. Only option for biLSTM models}
#'  }
#' @param geometric.theta angle to use for geometric transformation
#' 
#' @export
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom immApex propertyEncoder onehotEncoder geometricEncoder getIR
#' 
#' @return Ibex encoded values from the autoencoder
Ibex.matrix <- function(input.data, 
                        chain = "Heavy", 
                        method = "encoder",
                        encoder.model = "VAE", 
                        encoder.input = "AF",
                          geometric.theta = pi/3) {
    #Chain Check
    if(chains %!in% c("Heavy", "Light")) {
      stop("Please select one of the following chains: 'Heavy', 'Light'")
    }
  
    #TODO Pair these down to the final models offered
    #Model Check
    valid.encoder.inputs <- c("atchleyFactors", "crucianiProperties", "FASGAI", "kideraFactors", "MSWHIM", "ProtFP", "stScales", "tScales", "VHSE", "zScales", "OHE")
    if(encoder.input %!in% valid.encoder.inputs) {
      stop("Please select one of the valid encoder inputs.")
    }
  
    #Will be used to filter output of getIR()
    loci <- ifelse(chain == "Heavy", "IGH", c("IGK", "IGL"))
    expanded.sequences <- grepl(".Exp", encoder.model)
    
    #Set Dictionary for embedding
    if(expanded.sequences) {
      dictionary <- c(amino.acids, "_")
    } else {
      dictionary <- amino.acids
    }
    
    #Getting Sequences
    BCR <- getIR(input.data, chain, sequence.type = "aa")[[1]]
    
    #Filtering out NA values 
    if(any(is.na(BCR[,2]))) {
      BCR <- BCR[-which(is.na(BCR[,2])),]
    }
    #Filtering out sequences that do not match gene locid
    if(any(!grep(paste0(loci, collapse = "|"), BCR[,"v"]))) {
      BCR <- BCR[-!grep(paste0(loci, collapse = "|"), BCR[,"v"]),]
    }
    
    #Checking Sequences - needs to be < 45 or 90
    checkLength(x = BCR[,2], expanded = expanded.sequences)
    length.to.use <- ifelse(expanded.sequences, 90, 45)
    
    
    if (method == "encoder") {
      print("Encoding Sequences...")
      
      if(encoder.model == "biLSTM") {
        encoded.values <- suppressMessages(onehotEncoder(BCR[,2],
                                                         max.length = length.to.use,
                                                         convert.to.matrix = FALSE,
                                                         sequence.dictionary = dictionary,
                                                         padding.symbol = "."
                                                        ))
      } else {
        if(encoder.input == "OHE") {
          encoded.values <- suppressMessages(onehotEncoder(BCR[,2],
                                                           max.length = length.to.use,
                                                           convert.to.matrix = TRUE,
                                                           sequence.dictionary = dictionary,
                                                           padding.symbol = "."))
        } else {
          encoded.values <- suppressMessages(propertyEncoder(BCR[,2], 
                                                             max.length = length.to.use,
                                                             method.to.use = encoder.input,
                                                             convert.to.matrix = TRUE))
        }
      }
      
      print("Calculating Latent Dimensions...")
      #Getting Model
      aa.model <- aa.model.loader(chain, encoder.input, encoder.model)
      reduction <- stats::predict(aa.model, encoded.values, verbose = 0)
        
    } else if (method == "geometric") {
        print("Performing geometric transformation...")
        reduction <- suppressMessages(geometricEncoder(BCR, theta = geometric.theta))
    }
    #TODO Check the barcode usage here
    reduction <- as.data.frame(reduction)
    barcodes <- reduction[,1]
    reduction <- reduction[,-1]
    rownames(reduction) <- barcodes
    colnames(reduction) <- paste0("Ibex_", seq_len(ncol(reduction)))
    return(reduction)
}

#' Ibex single cell calculation
#'
#' Run Ibex algorithm with Seurat or SingleCellExperiment pipelines
#'
#' @examples
#' ibex_example <- runIbex(ibex_example, 
#'                         chains = "Heavy",
#'                         method = "encoder",
#'                         encoder.model = "VAE",
#'                         encoder.input = "AF")
#'                         
#' @param sc.data Single Cell Object in Seurat or SingleCell Experiment format
#' @param chain Heavy or Light
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on the BLOSUM62 matrix
#' @param encoder.model "CNN" = convolutional neural network-based autoencoder, 
#' "VAE" = variation autoencoder, or "biLSTM" for Bidirectional LSTM
#' @param encoder.input Type of model to use
#'   \itemize{
#'     \item{Amino Acid Properties: "atchleyFactors", "crucianiProperties", "FASGAI", "kideraFactors", "MSWHIM", "ProtFP", "stScales", "tScales", "VHSE", "zScales"}
#'     \item{"OHE" for One Hot Encoding. Only option for biLSTM models}
#'  }
#' @param geometric.theta angle to use for geometric transformation
#' @param reduction.name Keyword to save Ibex reduction. Useful if you want
#' to try Ibex with multiple parameters 
#' @export
#' @return Seurat or SingleCellExperiment object with Ibex dimensions placed 
#' into the dimensional reduction slot. 
#' 
runIbex <- function(sc.data, 
                    chains = "Heavy", 
                    method = "encoder",
                    encoder.model = "VAE", 
                    encoder.input = "AF",
                    geometric.theta = pi,
                    reduction.name = "Ibex") {
    checkSingleObject(sc.data)
    #TODO add check for CTAA columns
    cells.chains <- rownames(sc.data[[]][!is.na(sc.data[["CTaa"]]),])
    sc.data <- subset(sc.data, cells = cells.chains)
    reduction <- Ibex.matrix(input.data = sc.data,
                             chains = chains, 
                             method = method,
                             encoder.model = encoder.model, 
                             encoder.input = encoder.input,
                             geometric.theta = geometric.theta)
    BCR <- getIR(input.data, chain, sequence.type = "aa")[[1]]
    sc.data <- adding.DR(sc.data, reduction, reduction.name)
    return(sc.data)
}
