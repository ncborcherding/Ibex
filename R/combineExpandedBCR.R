#' combineBCR for CDR1/2/3 sequences
#'
#' @details This function enhances BCR processing by incorporating additional 
#' sequence information from CDR1 and CDR2 regions before applying the BCR 
#' combination logic. The function depends on 
#' \code{\link[scRepertoire]{combineBCR}} from the scRepertoire package.
#'
#' @param input.data List of filtered contig annotations.
#' @param samples Character vector. Labels of samples (required).
#' @param ID Character vector. Additional sample labeling (optional).
#' @param call.related.clones Logical. Whether to call related clones based on 
#'   nucleotide sequence and V gene. Default is `TRUE`.
#' @param threshold Numeric. Normalized edit distance for clone clustering. 
#' Default is `0.85`.
#' @param removeNA Logical. Whether to remove any chain without values. Default 
#' is `FALSE`.
#' @param removeMulti Logical. Whether to remove barcodes with more than two 
#' chains. Default is `FALSE`.
#' @param filterMulti Logical. Whether to select the highest-expressing light 
#' and heavy chains. Default is `TRUE`.
#' @param filterNonproductive Logical. Whether to remove nonproductive chains. 
#' Default is `TRUE`.
#'
#'@return A list of consolidated BCR clones with expanded CDR sequences.
#' @seealso 
#' \code{\link[scRepertoire]{combineBCR}}
#'
#' @importFrom scRepertoire combineBCR
#' @export
combineExpandedBCR <- function(input.data,
                               samples = NULL,
                               ID = NULL,
                               call.related.clones = TRUE,
                               threshold = 0.85,
                               removeNA = FALSE,
                               removeMulti = FALSE,
                               filterMulti = TRUE,
                               filterNonproductive = TRUE) {
  
  # Ensure input is a list of data frames
  if (!is.list(input.data) || !all(sapply(input.data, is.data.frame))) {
    stop("Input data must be a list of data frames.")
  }
  
  # Modify each data frame in the list
  modified_data <- lapply(input.data, function(df) {
    if (!all(c("cdr1", "cdr2", "cdr3") %in% colnames(df))) {
      stop("Each data frame must contain 'cdr1', 'cdr2', and 'cdr3' columns.")
    }
    
    # Create concatenated CDR sequence
    df$cdr3 <- paste(df$cdr1, df$cdr2, df$cdr3, sep = "-")
    df$cdr3_nt<- paste(df$cdr1_nt, df$cdr2_nt, df$cdr3_nt, sep = "-")
    
    return(df)
  })
  
  # Call combineBCR() on the modified data
  combined_result <- combineBCR(input.data = modified_data,
                                samples = samples,
                                ID = ID,
                                call.related.clones = call.related.clones,
                                threshold = threshold,
                                removeNA = removeNA,
                                removeMulti = removeMulti,
                                filterMulti = filterMulti,
                                filterNonproductive = filterNonproductive)
  
  return(combined_result)
}
