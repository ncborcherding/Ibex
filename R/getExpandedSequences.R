#' Extract full CDR sequences from contig files
#' 
#' This function is to generate CDR1/2/3 expanded
#' sequences that are functional with the appropriate model.
#' 
#' #' @examples
#' x <- quietBCRgenes(ibex_example)
#' 
#' @param sc A list of contig files with the "cdr1", "cdr2" and "cdr3" columns
#' @param chain Heavy or Light chain to extract
#' @importFrom dplyr %>% group_by slice_max
#' @export
#' @return A reorganized and filtered contig list with expanded_aa sequences


getExpandedSequences <- function(contig.files,
                                 chain = "Heavy") {
  
  if(chain == "Heavy") {
    loci <- "IGH"
  } else if (chain == "Light") {
    loci <- c("IGK, IGL")
  }
  
  #Check for cdr regions
  if(!all(c("cdr1", "cdr2", "cdr3") %in% colnames(contig.files[[1]]))) {
    stop("cdr1, cdr2, and/or cdr3 is missing from the contig files")
  }
  
  lapply(contig.files, function(x) {
    x <- x[!is.na(x[,"cdr3"]),]
    x <- x[x$chain %in% loci, ]
    x <- x %>%
      group_by(barcode) %>%
      slice_max(order_by = reads, n = 1) %>%
      as.data.frame()
    
    tmp <- data.frame(barcode = x$barcode, 
                      expanded_aa = paste0(x$cdr1, "_", x$cdr2, "_", x$cdr3), 
                      v = x$v_gene,
                      d = x$d_gene, 
                      j = x$j_gene,
                      c = x$c_gene)
    
  }) -> contig.files
  return(contig.files)
}
