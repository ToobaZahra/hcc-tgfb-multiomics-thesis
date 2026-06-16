# 15_methylation_analysis.R
library(SummarizedExperiment)
library(limma)

# Load methylation data
meth_se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")

cat("Dimensions:", dim(meth_se), "\n")
cat("Sample types:\n")
print(table(meth_se$shortLetterCode))

# Keep only tumor (TP) and normal (NT)
keep <- meth_se$shortLetterCode %in% c("TP", "NT")
meth_se <- meth_se[, keep]
cat("After filtering - Tumor:", sum(meth_se$shortLetterCode == "TP"),
    "Normal:", sum(meth_se$shortLetterCode == "NT"), "\n")

# Extract beta matrix
beta_mat <- assay(meth_se)
cat("Beta matrix dimensions:", dim(beta_mat), "\n")

# Remove probes with >20% missing values
na_frac <- rowMeans(is.na(beta_mat))
beta_mat <- beta_mat[na_frac < 0.2, ]
cat("After NA filter:", nrow(beta_mat), "probes remaining\n")

# Remove probes with low variance
probe_var <- rowVars(beta_mat, na.rm = TRUE)
beta_mat <- beta_mat[probe_var > 0.01, ]
cat("After variance filter:", nrow(beta_mat), "probes remaining\n")

# Save filtered matrix
saveRDS(beta_mat, "data/raw/tcga_lihc_meth_filtered.rds")
cat("Filtered beta matrix saved.\n")

# Differential methylation with limma
sample_type <- factor(meth_se$shortLetterCode, levels = c("NT", "TP"))
design <- model.matrix(~ sample_type)

# Impute missing values with row mean
beta_imputed <- t(apply(beta_mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

fit <- lmFit(beta_imputed, design)
fit <- eBayes(fit)
results_meth <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

cat("Differential methylation done.\n")
cat("Hypermethylated (delta-beta > 0):", sum(results_meth$logFC > 0.1 & results_meth$adj.P.Val < 0.05), "\n")
cat("Hypomethylated (delta-beta < 0):", sum(results_meth$logFC < -0.1 & results_meth$adj.P.Val < 0.05), "\n")

saveRDS(results_meth, "03_results/dmeth_tcga_tumor_vs_normal.rds")
write.csv(results_meth, "03_results/dmeth_tcga_tumor_vs_normal.csv", row.names = TRUE)
cat("Results saved.\n")