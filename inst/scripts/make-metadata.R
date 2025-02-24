#!/usr/bin/env Rscript
## make-metadata.R
## This script scans inst/extdata, builds a metadata.csv file
## for ExperimentHub or AnnotationHub submission.


PKG_NAME <- "Ibex"       
BIOC_VERSION <- "3.21"            
MAINTAINER <- "Nick Borcherding <ncborch@gail.com.com>"
DATA_PROVIDER <- "Consolidated Sources: IReceptor, OAS, and GEO"
SOURCE_URL <-  "https://github.com/BorchLab"
SOURCE_VERSION <- NA         
SOURCE_TYPE <- "CSV"              
GENOME <- NA                      
COORDINATE_1_BASED <- NA         
DESCRIPTION <- "Keras-based deep learning encoder for BCR sequences."


# 2) Locate the data files in inst/extdata
path_to_extdata <- file.path("inst", "extdata")
files <- list.files(path_to_extdata, full.names = TRUE, pattern = ".keras")

# 3) Helper function: guess DispatchClass and RDataClass from file extension
inferDispatchClass <- function(file_ext) {
  switch(
    tolower(file_ext),
    "rds" = "Rds",
    "rda" = "Rda",
    "csv" = "FilePath",
    "tsv" = "FilePath",
    "txt" = "FilePath",
    "FilePath"
  )
}

inferRDataClass <- function(dispatchClass) {
  # Adjust to reflect how your data is actually loaded in R.
  if (dispatchClass %in% c("Rds", "Rda")) {
    return("SummarizedExperiment")  # or whatever class your objects are
  } else {
    return("character") # or NA, if you just return a path
  }
}

# 4) Build metadata data.frame row-by-row
metadata_list <- lapply(files, function(f) {
  # Example: f == "inst/extdata/somefile.rds"
  file_name <- basename(f)
  file_ext <- tolower(tools::file_ext(f))  # "rds", "rda", "csv", etc.
  
  dispatchClass <- inferDispatchClass(file_ext)
  rDataClass   <- inferRDataClass(dispatchClass)
  
  # The Title could simply be the file name or something more descriptive
  title <- file_name  
  components <- stringr::str_split(title, "_")[[1]]
  
  # Adaptive Variables
  description <- paste0(DESCRIPTION, 
                        " Chain: ", components[2], 
                        ", Architecture: ", components[3], 
                        ", Encoding Method: ", components[4])
  SPECIES <- ifelse(grepl("Human", title), "Homo sapiens", "Mus musculus")
  TAXONOMY_ID <- ifelse(grepl("Human", title), "9606", "10090")
  rDataPath <- paste0("records/14919286/files/", file_name)
  # We assemble a named vector or list for each file:
  c(
    Title = title,
    Description = description,
    BiocVersion = BIOC_VERSION,
    Genome = as.character(GENOME),
    SourceType = SOURCE_TYPE,
    SourceUrl = SOURCE_URL,
    SourceVersion = SOURCE_VERSION,
    Species = SPECIES,
    TaxonomyId = TAXONOMY_ID,
    Coordinate_1_based = ifelse(is.na(COORDINATE_1_BASED), NA, 
                                as.character(COORDINATE_1_BASED)),
    DataProvider = DATA_PROVIDER,
    Maintainer = MAINTAINER,
    RDataClass = rDataClass,
    DispatchClass = dispatchClass,
    Location_Prefix = "https://zenodo.org/",
    RDataPath = rDataPath,
    Tags = paste("BCR", "scRNA-seq", "Encoder", "Model", sep = ":")
  )
})

# 5) Convert this list of named vectors to a data.frame
metadata_df <- do.call(rbind, lapply(metadata_list, as.data.frame.list))

# 6) Write out the metadata.csv to inst/extdata
output_csv <- file.path("inst", "extdata", "metadata.csv")
write.csv(metadata_df, file = output_csv, row.names = FALSE, quote = TRUE)

