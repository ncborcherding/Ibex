#' Cluster clones using the Ibex dimensional reductions
#' 
#' Use this to return clusters for clonotypes based on 
#' the \link[bluster]{bluster} clustering parameters.
#' 
#' @examples
#' \dontrun{
#' sc <- clonalCommunity(sc, 
#'                       reduction.name = NULL, 
#'                       cluster.parameter = KNNGraphParam())
#' }
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format. In addition, the outputs of distReduction()
#' and aaReduction() can be used.
#' @param reduction.name Name of the dimensional reduction output from runIbex()
#' @param cluster.parameter The community detection algorithm in \link[bluster]{bluster}
#' @param ... For the generic, further arguments to pass to specific methods.
#' @importFrom bluster clusterRows NNGraphParam HclustParam KmeansParam KNNGraphParam PamParam SNNGraphParam SomParam DbscanParam
#' @importFrom igraph simplify spectrum graph_from_edgelist E `E<-`
#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @return Single-Cell Object with ibex.clusters in the meta.data
clonalCommunity <- function(sc, 
                            reduction.name = NULL, 
                            cluster.parameter=KNNGraphParam(k=30, ...), 
                            ...) {
  if (inherits(x=sc, what ="Seurat")) { 
    dim.red <- sc[[reduction.name]] 
    dim.red <- dim.red@cell.embeddings
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
    dim.red <- reducedDim(sc, reduction.name)
  } else {
    if(inherits(x=sc, what ="dist")) {
      mat <- sc
      #mat[is.na(mat)] <- 0
      dimension <- attr(mat, "Size")
      edge <- NULL
      for (j in seq_len(dimension)[-1]) {
        row <- dist.convert(mat,j)
        tmp.edge <- data.frame("from" = j, "to" = seq_len(j)[-j], weight = row)
        edge <- rbind(edge, tmp.edge)
      }
      edge <- na.omit(edge)
      g <- graph.edgelist(as.matrix(edge[,c(1,2)]), directed = FALSE)
      E(g)$weights <- edge$weight
      g <- simplify(g)
      eigen <- spectrum(g, 
                        which = list(howmany = 30), 
                        algorithm = "arpack")
      dim.red <- eigen$vectors
    } else {
      dim.red <- sc
    }
  }
  clusters <- suppressWarnings(clusterRows(dim.red, BLUSPARAM=cluster.parameter))
  clus.df <- data.frame("ibex.clusters" = paste0("ibex.", clusters))
  if (inherits(x=sc, what ="Seurat") | inherits(x=sc, what ="SingleCellExperiment")) {
    rownames(clus.df) <- rownames(dim.red)
    sc <- add.meta.data(sc, clus.df, colnames(clus.df))
    return(sc)
  } 
  return(clus.df)
  
}