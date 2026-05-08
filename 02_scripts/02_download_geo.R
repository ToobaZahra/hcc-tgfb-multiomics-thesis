# 02_download_geo.R
# Download GEO replication datasets via GEOquery

library(GEOquery)

geo_dir <- "data/raw/geo"
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

geo_ids <- c("GSE14520", "GSE36376", "GSE76427")

for (gse in geo_ids) {
  cat("\n=== Downloading", gse, "===\n")
  gse_obj <- getGEO(gse, GSEMatrix = TRUE, destdir = geo_dir, AnnotGPL = TRUE)
  saveRDS(gse_obj, file = file.path(geo_dir, paste0(gse, ".rds")))
  cat("Saved", gse, "with", length(gse_obj), "platform(s)\n")
}

# Verify after download
cat("\n=== Verification ===\n")
for (gse in geo_ids) {
  obj <- readRDS(file.path(geo_dir, paste0(gse, ".rds")))
  for (i in seq_along(obj)) {
    cat(gse, "platform", i, "samples:", ncol(obj[[i]]), "\n")
  }
}