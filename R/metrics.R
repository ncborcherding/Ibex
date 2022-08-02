#Transforming the CDR3 sequence into vectors
#' @importFrom keras load_model_hdf5
#' @importFrom reticulate array_reshape
aaProperty <- function(TCR, 
                       AA.properties = AA.properties) { 
  return <- list() ### Need to add reference data
  reference <- Trex.Data[[1]] #AA properties
  col.ref <- grep(tolower(paste(AA.properties, collapse = "|")), colnames(reference))
  length <- NULL
  if (AA.properties == " both") {
    column.ref <- unique(sort(c(AF.col, KF.col)))
  } else {
    column.ref <- unique(sort(col.ref))
  }
  chain <- names(TCR)
  for (i in seq_along(TCR)) {
    membership <- TCR[[i]]
    names <- membership$barcode
    array.reshape <- NULL
    aa.model <- quiet(aa.model.loader(chain[[i]], AA.properties))
    range <- aa.range.loader(chain[[i]], AA.properties, Trex.Data) 
    local.min <- range[[1]]
    local.max <- range[[2]]
    }
    cells <- unique(membership[,"barcode"])
    score <- NULL
    for (n in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[n],]$cdr3_aa
      refer <- unlist(strsplit(tmp.CDR, ""))
      refer <- c(refer, rep(NA, 60 - length(refer)))
      if(AA.properties == "OHE") {
        int <- one.hot.organizer(refer)
        array.reshape.tmp <- array_reshape(int, 1260)
      }else {
        int <- reference[match(refer, reference$aa),c(1,col.ref)]
        int <- as.matrix(int[,-1])
        array.reshape.tmp <- array_reshape(int, length(col.ref)*60)
      }
      score.tmp <- auto.embedder(array.reshape.tmp, aa.model, local.max, local.min, AA.properties)
      score <- rbind(score, score.tmp)
    }
    score <- data.frame(unique(membership[,"barcode"]), score)
    barcodes <- score[,1]
    score <- score[,-1]
    rownames(score) <- barcodes
    colnames(score) <- paste0("Trex_", seq_len(ncol(score)))
    return(score)
}



