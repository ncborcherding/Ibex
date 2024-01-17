getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = "")) # could move testdata 1 dir lvl up nstead
}