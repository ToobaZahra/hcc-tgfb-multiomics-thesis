# 14_download_meth.R
library(TCGAbiolinks)

# Query TCGA-LIHC methylation 450K
query_meth <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

# Download
GDCdownload(query_meth, method = "api", files.per.chunk = 10,
            directory = "data/raw/GDCdata")

# Prepare
meth_se <- GDCprepare(query_meth, directory = "data/raw/GDCdata")
saveRDS(meth_se, "data/raw/tcga_lihc_meth450_se.rds")
cat("Methylation data saved.\n")
cat("Dimensions:", dim(meth_se), "\n")