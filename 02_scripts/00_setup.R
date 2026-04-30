# 00_setup.R - Run once to install all required packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.20", ask = FALSE)

bioc_pkgs <- c(
  "TCGAbiolinks", "DESeq2", "limma", "edgeR",
  "minfi", "ChAMP", "clusterProfiler", "org.Hs.eg.db",
  "GEOquery", "EnhancedVolcano", "ComplexHeatmap"
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

cran_pkgs <- c("tidyverse", "survival", "survminer", "glmnet", "shiny")
install.packages(cran_pkgs)

print(sessionInfo())