# 01_download_tcga.R
# Download TCGA-LIHC RNA-seq (tumor + normal).
# Output goes to data/raw/ which is gitignored.

library(TCGAbiolinks)
library(SummarizedExperiment)

data_dir <- "data/raw"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# RNA-seq: STAR counts for tumor + normal
rna_query <- GDCquery(
  project       = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(rna_query, directory = data_dir, method = "api", files.per.chunk = 20)
rna_se <- GDCprepare(rna_query, directory = data_dir)
saveRDS(rna_se, file = file.path(data_dir, "tcga_lihc_rnaseq_se.rds"))

cat("RNA-seq SummarizedExperiment saved.\n")
cat("Dimensions:", dim(rna_se), "\n")
cat("Sample types:\n")
print(table(rna_se$sample_type))

# Clinical data (small, fast)
clin <- GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")
saveRDS(clin, file = file.path(data_dir, "tcga_lihc_clinical.rds"))
cat("Clinical data rows:", nrow(clin), "\n")