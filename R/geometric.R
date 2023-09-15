.blosum62 <- matrix(
  c(
    4, -2, 0, -2, -1, -2, 0, -2, -1, -1, -1, -2, -1, -1, -1, 1, 0, -3, -2, 0,
    -2, 4, -3, 4, 1, -3, -1, 0, -3, 0, 0, -3, 3, 0, -3, 0, -1, -3, -2, -1,
    0, -3, 9, -3, -4, -2, -3, -3, -1, -3, -3, 1, 0, -3, 0, 1, 0, -2, -2, -3,
    -2, 4, -3, 6, 2, -3, -1, -1, -3, -1, -1, -4, 1, -1, -3, 0, -1, -4, -3, -1,
    -1, 1, -4, 2, 5, -3, -2, 0, -3, 1, 2, -3, 4, 1, -3, -1, -1, -3, -2, -1,
    -2, -3, -2, -3, -3, 6, -3, -1, 0, -3, -4, 0, -3, -3, -1, 0, -2, 1, 3, -3,
    0, -1, -3, -1, -2, -3, 6, -2, -4, -2, -3, 0, -2, -2, -2, 0, -2, -3, -2, -1,
    -2, 0, -3, -1, 0, -1, -2, 8, -3, -1, -1, -2, 1, 0, -2, 0, -2, -3, -2, 0,
    -1, -3, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -1, -1, -1, -2, -3, -1, 3,
    -1, 0, -3, -1, 1, -3, -2, -1, -3, 5, 1, -3, 1, 3, -3, -2, -1, -3, -1, 1,
    -1, 0, -3, -1, 2, -4, -3, -1, 2, 1, 5, -2, 0, 1, -1, -1, -1, -2, -1, 1,
    -2, -3, 1, -4, -3, 0, 0, -2, 1, -3, -2, 6, -3, -4, 0, 1, 0, -4, -3, -2,
    -1, 3, 0, 1, 4, -3, -2, 1, -3, 1, 0, -3, 6, 0, -3, 0, -1, -3, -2, -1,
    -1, 0, -3, -1, 1, -3, -2, 0, -1, 3, 1, -4, 0, 5, -2, -1, -1, -2, -1, 1,
    -1, -3, 0, -3, -3, -1, -2, -2, -1, -3, -1, 0, -3, -2, 7, -1, -2, -4, -3, -2,
    1, 0, 1, 0, -1, 0, 0, 0, -1, -2, -1, 1, 0, -1, -1, 4, 1, -3, -2, 0,
    0, -1, 0, -1, -1, -2, -2, -2, -2, -1, -1, 0, -1, -1, -2, 1, 5, -2, -2, 0,
    -3, -3, -2, -4, -3, 1, -3, -3, -3, -3, -2, -4, -3, -2, -4, -3, -2, 11, 2, -3,
    -2, -2, -2, -3, -2, 3, -2, -2, -1, -1, -1, -3, -2, -1, -3, -2, -2, 2, 7, -1,
    0, -1, -3, -1, -1, -3, -1, 0, 3, 1, 1, -2, -1, 1, -2, 0, 0, -3, -1, 4
  ), nrow = 20, byrow = TRUE,
  dimnames = list(c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
                  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"))
)

.aa_to_blosum62 <- function(aa_sequence) {
  aa_map <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  aa_indices <- match(strsplit(as.character(aa_sequence), '')[[1]], aa_map)
  return(.blosum62[aa_indices, ])
}

# Main function to encode each sequence
#' @importFrom matrixStats colWeightedMeans
.encode_sequence <- function(sequence, dim, rotation_matrix) {
  blosum_vectors <- .aa_to_blosum62(sequence)
  
  # Initialize an empty matrix to store the transformed points
  transformed_points <- matrix(nrow=0, ncol=dim)
  
  # Apply the unitary transformation
  for (i in 1:(dim(blosum_vectors)[1])) {
    point <- t(blosum_vectors[i, ])
    
    # Apply the rotation matrix
    transformed_point <- rotation_matrix %*% matrix(point, nrow=dim, ncol=1)
    
    # Collect the transformed points
    transformed_points <- rbind(transformed_points, t(transformed_point))
  }
  
  # Average the transformed points to get a single 20D point
  avg_point <- colWeightedMeans(transformed_points)
  
  return(avg_point)
}

.geometric.encoding <- function(BCR, 
                                theta =theta) {
  dim <- 20
  rotation_matrix <- diag(dim)
  
  # Generate ten 2D rotation matrices and insert them into the 20D rotation matrix
  rotation_2d = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol=2)
  
  for (i in seq(1, dim, by=2)) {
    rotation_matrix[i:(i+1), i:(i+1)] = rotation_2d
  }
  
  transformed <- lapply(BCR[[1]]$cdr3_aa, function(x) {
    tmp <- .encode_sequence(x, dim, rotation_matrix)
    tmp
  })
  score <- do.call(rbind,transformed)
  rownames(score) <- BCR[[1]]$barcode
  colnames(score) <- paste0("Ibex_", seq_len(ncol(score)))
  return(score)
}


