#' Reduce a Single-Cell Object to Representative Cells
#'
#' This function generates a single-cell object with a reduced representation 
#' of RNA expression by clone. The approach is inspired by the method introduced 
#' in [CoNGA](https://pubmed.ncbi.nlm.nih.gov/34426704/). Users can 
#' generate either a mean representation of features by clone or identify a 
#' representative cell using count-based minimal Euclidean distance. 
#' Please read and cite the original work by the authors of CoNGA.
#'
#' @examples
#' #' # Get Data
#' ibex_example <- get(data("ibex_example"))
#' 
#' ibex.clones <- CoNGAfy(ibex_example, 
#'                        method = "dist")
#'
#' ibex.clones <- CoNGAfy(ibex_example, 
#'                        method = "mean")
#'
#' @param input.data A single-cell dataset in Seurat or SingleCellExperiment format.
#' @param method Character. Specifies the method to reduce the dataset:
#'   - "mean" - Computes the mean expression of selected features across cells in each clonotype.
#'   - "dist" - Uses PCA reduction to identify the cell with the minimal Euclidean distance within each clonotype group.
#' @param features Character vector. Selected genes for the reduction. If `NULL` (default), all genes are used.
#' @param assay Character. The name of the assay or assays to include in the output. Defaults to the active assay.
#' @param meta.carry Character vector. Metadata variables to carry over from the input single-cell object to the output.
#'
#' @return A reduced single-cell object where each clonotype is represented by a single cell.
#'
#' @export
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<-
#' @importFrom SummarizedExperiment assay assay<- SummarizedExperiment colData<-
#'             colData

CoNGAfy <- function(input.data, 
                    method = "dist", 
                    features = NULL, 
                    assay = "RNA", 
                    meta.carry = c("CTaa", "CTgene")) {
    if(inherits(input.data, "Seurat")) {
      cells.chains <- rownames(input.data[[]][!is.na(input.data[["CTaa"]]),])
      input.data <- subset(input.data, cells = cells.chains)
    } else if (inherits(input.data, "SingleCellExperiment")) {
      cells.chains <- rownames(as.data.frame(colData(input.data)[!is.na(input.data$CTaa),]))
      input.data <- input.data[,which(colnames(input.data) %in% cells.chains)]
    } else {
      stop("The input.data is not a Seurat or SingleCellExperiment object.")
    }
    conga <- NULL
    if(method == "mean") {
        for (x in seq_along(assay)) {
            conga[[x]] <- .CoNGA.mean(input.data, features, assay[x])
            
        }
    } else if(method == "dist") {
        for (x in seq_along(assay)) {
            conga[[x]] <- .CoNGA.dist(input.data, features, assay[x])
            
        }
        
    }
    names(conga) <- assay
    if (inherits(x=input.data, what ="Seurat")) {
        sc.output <- CreateSeuratObject(conga[[1]], assay = names(conga)[1], project = "Ibex")
        if(length(conga) > 1) {
            for(y in 2:length(conga)) {
                sc.output[[names(conga)[y]]] <- CreateAssayObject(conga[[y]])
            }
        }
        CTge <- unique(input.data[[]][,c(meta.carry)])
    } else if (inherits(x=input.data, what ="SingleCellExperiment")) {
        sc.output <- SingleCellExperiment(assay = conga[[1]])
        if(length(conga) > 1) {
            for(y in 2:length(conga)) {
                altExp(sc.output, "BEAM") <- SummarizedExperiment(
                  assays = list(
                    counts = as.matrix(conga[[y]])
                  ),
                  colData = colData(sc.output)
                )
            }
        }
        sc.output$CTaa <- rownames(sc.output@colData)
        CTge <- data.frame(unique(input.data@colData[,c(meta.carry)]))
    }
    CTge <- CTge[!duplicated(CTge$CTaa),]
    clones <- unique(CTge$CTaa)
    rownames(CTge) <- clones
    colnames(CTge) <- c("CTaa", "CTgene")
    sc.output <- add.meta.data(sc.output, CTge, colnames(CTge))
    return(sc.output)
}

# Pulls Assay Data
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom SingleCellExperiment altExp
grabAssay <- function(input.data, assay) {
  if (inherits(x=input.data, what ="Seurat")) {
    data.use <- input.data[[assay]]$counts
  } else if (inherits(x=input.data, what ="SingleCellExperiment")){
    if(assay %in% assayNames(input.data) | assay == "RNA") {
      if(assay == "RNA")  assay <- "counts"
      data.use <- assay(input.data, name = assay)
      } else {
        data.use <- assay(altExp(input.data), name = assay)
      }
    }
  return(data.use)
}

# Calculate best representation individual clones
#' @importFrom SummarizedExperiment assay
#' @importFrom SeuratObject GetAssayData
#' @importFrom methods is
#' @importFrom stats dist
#' @keywords internal
#' @noRd
.CoNGA.dist <- function(input.data, 
                        features = NULL, 
                        assay = "RNA") {
  # Ensure 'assay' is character (vector or single string)
  if (!is.character(assay)) {
    stop("'assay' must be a character vector or a single character string.")
  }
  
  # Grab clone meta-information; here we assume 'grabMeta' returns a DataFrame or data.frame
  meta <- grabMeta(input.data)
  # Create a small table of CTaa assignments
  ct_col <- "CTaa"
  if (!ct_col %in% colnames(meta)) {
    stop("The metadata must contain a column named 'CTaa'.")
  }
  meta_ct <- data.frame(CTaa = meta[, ct_col], row.names = rownames(meta))
  
  # Identify number of cells per clone
  clone_tab <- table(meta_ct$CTaa)
  multi_clone_names <- names(clone_tab[clone_tab > 1])  # clones with >1 cell
  single_clone_names <- names(clone_tab[clone_tab == 1])# clones with exactly 1 cell
  
  # Function to process a single assay
  process_single_assay <- function(assay_name) {
    # Pull the correct data matrix from input.data
    data_mat <- grabAssay(input.data, assay_name)
    
    # Subset the features if requested
    features_to_use <- features %||% rownames(data_mat)
    features_to_use <- intersect(features_to_use, rownames(data_mat))
    
    # If no features remain, warn
    if (length(features_to_use) == 0) {
      warning("No overlapping features found in assay '", assay_name, "'. Returning empty matrix.")
      return(matrix(nrow = 0, ncol = 0))
    }
    
    # Subset 'data_mat' to only those features
    data_mat_use <- data_mat[features_to_use, , drop = FALSE]
    
    # We now find the "best representation" for each multi-cell clone by minimal sum of distances
    # Start with barcodes that are single-cell clones (they trivially represent themselves)
    best_barcodes <- rownames(meta_ct)[meta_ct$CTaa %in% single_clone_names]
    
    # For each multi-cell clone, compute distances and pick the cell with smallest total distance
    for (clone_name in multi_clone_names) {
      clone_cells <- rownames(meta_ct)[meta_ct$CTaa == clone_name]
      # Distances are among rows of data_mat_use
      dist_mat <- as.matrix(dist(t(as.matrix(data_mat_use[, clone_cells, drop = FALSE]))))
      
      # rowSums(dist_mat) is sum of distances from each cell to all others in the clone
      chosen_idx <- which.min(rowSums(dist_mat))
      chosen_cell <- clone_cells[chosen_idx]
      best_barcodes <- c(best_barcodes, chosen_cell)
    }
    
    # Finally, subset original matrix to these 'best_barcodes'
    data_return <- data_mat_use[, best_barcodes, drop = FALSE]
    # Rename columns to the clone name for clarity
    colnames(data_return) <- meta_ct$CTaa[match(best_barcodes, rownames(meta_ct))]
    
    return(data_return)
  }
  
  # If user passed multiple assays, return a list
  if (length(assay) > 1) {
    results_list <- lapply(assay, process_single_assay)
    names(results_list) <- assay
    return(results_list)
  } else {
    # If user passed a single assay, return a single matrix
    return(process_single_assay(assay))
  }
}

# Calculate mean across individual clones
#' @importFrom rlang %||%
#' @importFrom Matrix sparse.model.matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom SeuratObject GetAssayData
#' @importFrom stats as.formula
#' @keywords internal
#' @noRd
.CoNGA.mean <- function(input.data, 
                        features = NULL, 
                        assay = "RNA") {
  # Ensure 'assay' is character (vector or single string)
  if (!is.character(assay)) {
    stop("'assay' must be a character vector or a single character string.")
  }
  
  # Grab clone meta-information
  meta <- grabMeta(input.data)
  ct_col <- "CTaa"
  if (!ct_col %in% colnames(meta)) {
    stop("The metadata must contain a column named 'CTaa'.")
  }
  meta_ct <- data.frame(CTaa = meta[, ct_col], row.names = rownames(meta))
  
  # Remove rows with NA in CTaa
  meta_ct <- meta_ct[which(rowSums(is.na(meta_ct)) == 0), , drop = FALSE]
  # Convert CTaa to a factor
  meta_ct$CTaa <- as.factor(meta_ct$CTaa)
  
  # Construct a model matrix with no intercept
  # ~0 + CTaa means we get one column per level of CTaa
  category_matrix <- sparse.model.matrix(
    as.formula('~0+CTaa'), 
    data = meta_ct
  )
  
  # Precompute column sums and scale columns to sum to 1
  col_sums <- Matrix::colSums(category_matrix)
  # remove columns with zero count if any
  keep_cols <- which(col_sums > 0)
  category_matrix <- category_matrix[, keep_cols, drop = FALSE]
  col_sums <- col_sums[keep_cols]
  
  # scale columns so each column sums to 1
  category_matrix <- sweep(category_matrix, MARGIN = 2, STATS = col_sums, FUN = "/")
  
  # Function to process a single assay
  process_single_assay <- function(assay_name) {
    data_mat <- grabAssay(input.data, assay_name)
    
    # Subset features if requested
    features_to_use <- features %||% rownames(data_mat)
    features_to_use <- intersect(features_to_use, rownames(data_mat))
    
    if (length(features_to_use) == 0) {
      warning("No overlapping features found in assay '", assay_name, "'. Returning empty matrix.")
      return(matrix(nrow = 0, ncol = 0))
    }
    
    data_mat_use <- data_mat[features_to_use, , drop = FALSE]
    
    # Multiply by the category matrix to get mean expression per clone
    # For each feature, we do feature_values %*% category_matrix 
    # (since category_matrix has columns that are "per-clone" indicators).
    data_return <- data_mat_use %*% category_matrix
    
    # Rename columns to reflect the clone name(s)
    colnames(data_return) <- gsub("^CTaa", "", colnames(category_matrix))
    
    return(data_return)
  }
  
  # If multiple assays, return a list
  if (length(assay) > 1) {
    results_list <- lapply(assay, process_single_assay)
    names(results_list) <- assay
    return(results_list)
  } else {
    # Single assay: return just a single matrix
    return(process_single_assay(assay))
  }
}
