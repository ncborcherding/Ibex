#' Ibex Single-Cell Calculation
#'
#' This function applies the Ibex algorithm to single-cell data, integrating seamlessly with 
#' Seurat or SingleCellExperiment pipelines. The algorithm generates latent dimensions using 
#' deep learning or geometric transformations, storing the results in the dimensional reduction slot.
#'
#' @examples
#' # Using the encoder method with a variational autoencoder
#' ibex_example <- runIbex(ibex_example, 
#'                         chain = "Heavy",
#'                         method = "encoder",
#'                         encoder.model = "VAE",
#'                         encoder.input = "atchleyFactors")
#'
#' # Using the geometric method with a specified angle
#' ibex_example <- runIbex(ibex_example, 
#'                         chain = "Heavy",
#'                         method = "geometric",
#'                         geometric.theta = pi)
#'
#' @param sc.data A single-cell dataset, which can be:
#'   \itemize{
#'     \item A Seurat object
#'     \item A SingleCellExperiment object
#'   }
#' @param chain Character. Specifies the chain to analyze:
#'   \itemize{
#'     \item "Heavy" for the heavy chain
#'     \item "Light" for the light chain
#'   }
#' @param method Character. Algorithm to use for generating latent dimensions:
#'   \itemize{
#'     \item "encoder" - Uses deep learning autoencoders
#'     \item "geometric" - Uses geometric transformations based on the BLOSUM62 matrix
#'   }
#' @param encoder.model Character. Type of autoencoder model to use:
#'   \itemize{
#'     \item "CNN" - Convolutional Neural Network-based autoencoder
#'     \item "VAE" - Variational Autoencoder
#'   }
#' @param encoder.input Character. Input features for the encoder model:
#'   \itemize{
#'     \item Amino Acid Properties: "atchleyFactors", "crucianiProperties", "kideraFactors", "MSWHIM", "tScales"
#'     \item "OHE" - One Hot Encoding 
#'   }
#' @param geometric.theta Numeric. Angle (in radians) for geometric transformation. 
#'   Used only when \code{method = "geometric"}.
#' @param reduction.name Character. The name to assign to the dimensional reduction. 
#'   This is useful for running Ibex with multiple parameter settings and saving results 
#'   under different names.
#'
#' @return An updated Seurat or SingleCellExperiment object with Ibex dimensions added 
#' to the dimensional reduction slot.
#' @export
runIbex <- function(sc.data, 
                    chains = "Heavy", 
                    method = "encoder",
                    encoder.model = "VAE", 
                    encoder.input = "AF",
                    geometric.theta = pi,
                    reduction.name = "Ibex") {
    checkSingleObject(sc.data)
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
