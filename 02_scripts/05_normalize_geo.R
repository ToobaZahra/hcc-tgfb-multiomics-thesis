# 05_normalize_geo.R
library(GEOquery)
library(limma)

geo_dir <- "data/raw/geo"

# Load each GEO dataset
gse14520 <- readRDS(file.path(geo_dir, "GSE14520.rds"))
gse36376 <- readRDS(file.path(geo_dir, "GSE36376.rds"))
gse76427 <- readRDS(file.path(geo_dir, "GSE76427.rds"))

# Function to normalize — with optional log2 transform
normalize_geo <- function(gse_obj, gse_name, log_transform = FALSE) {
  eset <- gse_obj[[1]]
  mat <- exprs(eset)
  
  if (log_transform) {
    cat(gse_name, "— applying log2 transform\n")
    mat <- log2(mat - min(mat) + 1)
  }
  
  mat_norm <- normalizeBetweenArrays(mat, method = "quantile")
  cat(gse_name, "— dimensions:", nrow(mat_norm), "probes x", ncol(mat_norm), "samples\n")
  return(mat_norm)
}

mat14520 <- normalize_geo(gse14520, "GSE14520", log_transform = FALSE)
mat36376 <- normalize_geo(gse36376, "GSE36376", log_transform = FALSE)
mat76427 <- normalize_geo(gse76427, "GSE76427", log_transform = FALSE)

saveRDS(mat14520, file.path(geo_dir, "GSE14520_norm.rds"))
saveRDS(mat36376, file.path(geo_dir, "GSE36376_norm.rds"))
saveRDS(mat76427, file.path(geo_dir, "GSE76427_norm.rds"))

cat("All GEO datasets normalized and saved.\n")