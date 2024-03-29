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
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format or
#' the output of combineBCR() in scRepertoire
#' @param chains Heavy or Light
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on BLOSUM62 matrix
#' @param encoder.model "AE" = dense autoencoder or "VAE" = variation autoencoder
#' @param encoder.input "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or "OHE" for
#' One Hot Autoencoder
#' @param theta angle to use for geometric transformation
#' 
#' @export
#' @importFrom SeuratObject CreateDimReducObject
#' 
#' @return Ibex encoded values from the autoencoder
Ibex.matrix <- function(sc, 
                        chains = "Heavy", 
                        method = "encoder",
                        encoder.model = "VAE", 
                        encoder.input = "AF",
                        theta = pi) {
    chains <- chain.checker(chains)
    BCR <- getBCR(sc, chains)
    checkLength(BCR[[1]])
    if (method == "encoder" && encoder.input %in% c("AF", "KF", "both", "all", "OHE")) {
        print("Calculating the encoding values...")
        reduction <- .encoder(BCR, encoder.input, encoder.model)
    } else if (method == "geometric") {
        print("Performing geometric transformation...")
        reduction <- .geometric.encoding(BCR, theta)
    }
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
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains Heavy or Light
#' @param chains Heavy or Light
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on BLOSUM62 matrix
#' @param encoder.model "AE" = dense autoencoder or "VAE" = variation autoencoder
#' @param encoder.input "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or "OHE" for
#' One Hot Autoencoder
#' @param theta angle to use for geometric transformation
#' @param reduction.name Keyword to save Ibex reduction. Useful if you want
#' to try Ibex with multiple parameters 
#' @export
#' @return Seurat or SingleCellExperiment object with Ibex dimensions placed 
#' into the dimensional reduction slot. 
#' 
runIbex <- function(sc, 
                    chains = "Heavy", 
                    method = "encoder",
                    encoder.model = "VAE", 
                    encoder.input = "AF",
                    theta = pi,
                    reduction.name = "Ibex") {
    checkSingleObject(sc)
    cells.chains <- rownames(sc[[]][!is.na(sc[["CTaa"]]),])
    sc <- subset(sc, cells = cells.chains)
    reduction <- Ibex.matrix(sc = sc,
                             chains = chains, 
                             method = method,
                             encoder.model = encoder.model, 
                             encoder.input = encoder.input,
                             theta = theta)
    BCR <- getBCR(sc, chains)
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}
