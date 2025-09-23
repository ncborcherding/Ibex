library(scRepertoire)
library(Seurat)
library(dplyr)
library(SummarizedExperiment)
library(SingleCellExperiment)
##################################
#scRNA/ADT loading and processing
#################################

tmp <-  Read10X("~/data/filtered_feature_bc_matrix")

SeuratObj <- CreateSeuratObject(counts = tmp$`Gene Expression`)
beam_assay <- CreateAssayObject(counts = tmp$`Antigen Capture`)

SeuratObj[["BEAM"]] <- beam_assay
SeuratObj <- subset(SeuratObj, subset = nFeature_RNA > 100) 
SeuratObj  <- RenameCells(object = SeuratObj , new.names = paste0("BEAM.sample_", rownames(SeuratObj[[]])))
SeuratObj[["mito.genes"]] <- PercentageFeatureSet(SeuratObj, pattern = "^mt-")

#Filtering step
standev <- sd(log(SeuratObj$nFeature_RNA))*2.5 #cutting off above standard deviation of 2.5
mean <- mean(log(SeuratObj$nFeature_RNA))
cut <- round(exp(standev+mean))
SeuratObj <- subset(SeuratObj, subset = mito.genes < 10 & nFeature_RNA < cut)

#Processing and Adding Contig Info
contigs <- read.csv("~/data/2k_BEAM-Ab_Mouse_HEL_5pv2_2k_BEAM-Ab_Mouse_HEL_5pv2_vdj_b_filtered_contig_annotations.csv")
clones <- combineBCR(contigs, samples = "BEAM.sample", removeNA = TRUE)
SeuratObj <- combineExpression(clones, SeuratObj, cloneCall="aa")

#Subset only cells with BCR and Heavy Chain 
cell.idx <- intersect(which(!is.na(SeuratObj$CTaa)), which(!is.na(stringr::str_split(SeuratObj$CTaa, "_", simplify = TRUE)[,1])))
SeuratObj <- subset(SeuratObj, cells = colnames(SeuratObj)[cell.idx])

#Processing RNA
DefaultAssay(SeuratObj) <- 'RNA'
SeuratObj <- NormalizeData(SeuratObj, verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>% 
  quietBCRgenes() %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)

#Removing negative control + B Cells
DefaultAssay(SeuratObj) <- 'BEAM'
SeuratObj <- subset(SeuratObj, subset = `negative-control` < 100, slot = "counts")

#Processing BEAM
VariableFeatures(SeuratObj) <- rownames(SeuratObj[["BEAM"]])
SeuratObj <- NormalizeData(SeuratObj, 
                           normalization.method = 'CLR',
                           margin = 2, ) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE, reduction.name = 'apca')

DefaultAssay(SeuratObj) <- 'RNA' 
###################################
#Making Example Data Set for Ibex
#################################

# Subset nondominate clones + random sampling of dominant
set.seed(42)
cell.idx <- unique(c(which(!grepl("CANWDGDYW", SeuratObj$CTaa)), sample(seq_len(nrow(SeuratObj[[]])), 154)))

ibex_example <- SeuratObj
saveRDS(ibex_example, file = "Ibex_FullExample.rds")

# Forming Example Data set in SCE format
ibex_example <- subset(ibex_example, cells = colnames(ibex_example)[cell.idx])
PCA <- Embeddings(ibex_example[["pca"]])
APCA <- Embeddings(ibex_example[["apca"]])
BEAM_counts <- GetAssayData(ibex_example, slot = "counts", assay = "BEAM")[1:4,]
BEAM_data   <- GetAssayData(ibex_example, slot = "data",   assay = "BEAM")[1:4,]
ibex_example <- as.SingleCellExperiment(ibex_example)
altExp(ibex_example, "BEAM") <- SummarizedExperiment(
  assays = list(
    counts = as.matrix(BEAM_counts),
    logcounts = as.matrix(BEAM_data) 
  ),
  colData = colData(ibex_example)
)
reducedDim(ibex_example, "pca") <- PCA
reducedDim(ibex_example, "apca") <- APCA

#Saving the built-in data set
save(ibex_example, file = "ibex_example.rda", compress = "xz")
ibex_vdj <- contigs
save(ibex_vdj, file = "ibex_vdj.rda", compress = "xz")