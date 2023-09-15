# Transforming the CDR3 sequence into vectors
#' @importFrom keras load_model_hdf5
#' @importFrom reticulate array_reshape
.encoder <- function(BCR, 
                     encoder.input = encoder.input, 
                     encoder.model = encoder.model) { 
  return <- list() ### Need to add reference data
  reference <- ibex.data[[1]] #AA properties
  col.ref <- grep(tolower(paste(encoder.input, collapse = "|")), colnames(reference))
  length <- NULL
  if (encoder.input == " both") {
    column.ref <- unique(sort(c(AF.col, KF.col)))
  } else {
    column.ref <- unique(sort(col.ref))
  }
  chain <- names(BCR)
  aa.model <- quiet(aa.model.loader(chain, encoder.input, encoder.model))
  membership <- BCR[[1]]
  names <- membership$barcode
  array.reshape <- NULL
    
  cells <- unique(membership[,"barcode"])
  score <- NULL
  for (n in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[n],]$cdr3_aa
      refer <- unlist(strsplit(tmp.CDR, ""))
      refer <- c(refer, rep(NA, 70 - length(refer)))
      if(encoder.input == "OHE") {
        int <- one.hot.organizer(refer)
        array.reshape.tmp <- array_reshape(int, 1470)
      }else {
        int <- reference[match(refer, reference$aa),c(1,col.ref)]
        int <- as.matrix(int[,-1])
        array.reshape.tmp <- array_reshape(int, length(col.ref)*70)
      }
      score.tmp <- auto.embedder(array.reshape.tmp, aa.model, encoder.input)
      score <- rbind(score, score.tmp)
    }
    score <- data.frame(unique(membership[,"barcode"]), score)
    barcodes <- score[,1]
    score <- score[,-1]
    rownames(score) <- barcodes
    colnames(score) <- paste0("Ibex_", seq_len(ncol(score)))
    return(score)
}



