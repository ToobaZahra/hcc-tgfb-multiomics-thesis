# 03_build_dds.R
# Build DESeq2 dataset from TCGA SummarizedExperiment

library(DESeq2)
library(SummarizedExperiment)

rna_se <- readRDS("data/raw/tcga_lihc_rnaseq_se.rds")

# Use unstranded counts assay
counts_mat <- assay(rna_se, "unstranded")
col_data <- colData(rna_se)

# Convert sample_type to factor (Normal as reference)
col_data$sample_type <- factor(col_data$sample_type,
                               levels = c("Solid Tissue Normal", "Primary Tumor"))

# Build DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ sample_type
)

cat("DDS dimensions:", dim(dds), "\n")
cat("Sample type table:\n")
print(table(dds$sample_type))

# Filter: keep genes with rowSums(counts) >= 10
keep <- rowSums(counts(dds)) >= 10
dds_filt <- dds[keep, ]

cat("Genes before filter:", nrow(dds), "\n")
cat("Genes after filter:",  nrow(dds_filt), "\n")

# Save
saveRDS(dds, file = "data/raw/tcga_lihc_dds_raw.rds")
saveRDS(dds_filt, file = "data/raw/tcga_lihc_dds_filtered.rds")