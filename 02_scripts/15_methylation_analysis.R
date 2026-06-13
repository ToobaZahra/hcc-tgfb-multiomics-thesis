# 15_methylation_analysis.R
library(SummarizedExperiment)
library(limma)

# Load methylation data
meth_se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")

cat("Dimensions:", dim(meth_se), "\n")
cat("Sample types:\n")
print(table(meth_se$shortLetterCode))
