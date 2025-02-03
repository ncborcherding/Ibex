# Ibex
Using BCR sequences for graph embedding

[![R-CMD-check](https://github.com/BorchLab/Ibex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/Ibex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BorchLab/Ibex/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/Ibex?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.borch.dev/uploads/screpertoire/articles/ibex)

<img align="right" src="https://github.com/BorchLab/Ibex/blob/main/www/ibex_hex.png" width="352" height="352">

## Introduction
Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/BorchLab/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. However, utilization of clonal indices for more complex analyses are still lacking, specifically in using clonality in embedding of single-cells. To this end, we developed an R package that uses deep learning to vectorize BCR sequences using order or translating the sequence into amino acid properties. The sister package to this is [Trex](https://github.com/BorchLab/Trex) for embedding of TCR sequence data.

# System requirements 

Ibex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have had BCRs added using [scRepertoire](https://github.com/BorchLab/scRepertoire). Ibex has been tested on OS X and Windows platforms.

**keras** is necessary to use the autoencoder function (this includes the set up of the tensorflow environment in R):

```r
##Install keras
install.packages("keras")

##Setting up Tensor Flow
library(reticulate)
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(keras)
install_keras()
```

# Installation

To run Ibex, open R and install Ibex from github: 

```r
devtools::install_github("BorchLab/Ibex")
```

# Usage/Demos

Ibex should be able to be run in popular R-based single-cell workflows, including Seurat and Bioconductor/Single-Cell Experiment formats.

## Quick Start 

Check out this [vignette](https://www.borch.dev/uploads/screpertoire/articles/ibex) for a quick start tutorial. 

<img align="center" src="https://github.com/BorchLab/Ibex/blob/main/www/graphicalAbstract.png">



## Autoencoded Matrix

The Ibex algorithm allows users to select BCR-based metrics to return autoencoded values to be used in dimensional reduction. If single-cell objects are not filtered for B cells with BCR,  `Ibex.matrix()` will still return values, however IBEX_1 will be based on the disparity of BCR-containing and BCR-non-containing cells based on the Ibex algorithm. 

```r
library(Ibex)
my_ibex <- Ibex.matrix(singleObject)
```

## Seurat or Single-Cell Experiment

You can run Ibex within your Seurat or Single-Cell Experiemt workflow. **Importantly** `runIbex()` will automatically filter single-cells that do not contain BCR information in the meta data of the single-cell object. 

```r
seuratObj_Bonly <- runIbex(seuratObj, #The single cell object
                           chains = "Heavy", #Use of "Heavy" or "Light" 
                           AA.properties = c("AF", "KF", "both", "OHE"), 
                           reduction.name = "Ibex", #Name designation for 
                           #the vectors to be added to the single-cell object)
                   
seuratObj_Bonly <- runIbex(seuratObj, reduction.name = "Ibex")
```

## After Running Ibex

From here, you can generate a tSNE/UMAP using the Ibex values, similar to the PCA values based on variable gene expression.

```r
seuratObj <- RunTSNE(seuratObj, reduction = "Ibex",  reduction.key = "Ibex_")
seuratObj <- RunUMAP(seuratObj, reduction = "Ibex",  reduction.key = "Ibex_")
```

If using Seurat package, the Ibex embedding information and gene expression PCA can be used to find the [Weighted Nearest Neighbors](https://pubmed.ncbi.nlm.nih.gov/34062119/). Before applying the WNN approach, best practice would be to remove the BCR-related genes from the list of variable genes and rerunning the PCA analysis. 

### Recalculate PCA without BCR genes with quietBCRgenes() function in Ibex.
```r
seuratObj <- quietBCRgenes(seuratObj)
seuratObj <- RunPCA(seuratObj)
```

### Running WNN approach
```r
seuratObj <- FindMultiModalNeighbors(seuratObj, 
                                     reduction.list = list("pca", "Ibex"), 
                                     dims.list = list(1:30, 1:20), 
                                     modality.weight.name = "RNA.weight")
                                     
seuratObj <- RunUMAP(seuratObj, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")
```
## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/Ibex/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the internal **ibex_example** would 
be extremely helpful.

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/Ibex/issues).

#### [Pull Requests](https://github.com/BorchLab/Ibex/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation
More information on Ibex is available at our [Biorxiv preprint](https://www.biorxiv.org/content/10.1101/2022.11.09.515787v2). 
