# 28_meth_expr_correlation.R
library(SummarizedExperiment)
library(dplyr)

# Load data
meth_se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds")
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

# Check dimensions
cat("Methylation samples:", ncol(meth_se), "\n")
cat("Expression samples:", ncol(vsd), "\n")

# Match tumor samples
meth_barcodes <- substr(colnames(meth_se), 1, 12)
expr_barcodes <- substr(colnames(vsd), 1, 12)
common <- intersect(meth_barcodes, expr_barcodes)
cat("Common samples:", length(common), "\n")



# Subset to common samples
meth_keep <- match(common, meth_barcodes)
expr_keep <- match(common, expr_barcodes)

meth_sub <- meth_se[, meth_keep]
expr_sub <- vsd[, expr_keep]

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# Get probe annotation
probe_annot <- as.data.frame(rowData(meth_se))
probe_annot$probe_id <- rownames(probe_annot)
nexus_probes <- probe_annot[probe_annot$gene %in% nexus_genes, c("probe_id", "gene")]

# Beta matrix for nexus probes
beta_mat <- assay(meth_sub)[nexus_probes$probe_id, ]
expr_mat <- assay(expr_sub)

# For each gene — correlate mean methylation vs expression
corr_results <- data.frame()

for (gene in nexus_genes) {
  # Get probes for this gene
  probes <- nexus_probes$probe_id[nexus_probes$gene == gene]
  probes <- probes[probes %in% rownames(beta_mat)]
  if (length(probes) == 0) next
  
  # Mean methylation across probes
  mean_meth <- colMeans(beta_mat[probes, , drop = FALSE], na.rm = TRUE)
  
  # Get expression
  ens <- res_df$ensembl_id[res_df$symbol == gene & !is.na(res_df$symbol)]
  if (length(ens) == 0) next
  ens <- ens[1]
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) == 0) next
  expr_vals <- expr_mat[row_match[1], ]
  
  # Correlation
  ct <- cor.test(mean_meth, expr_vals, method = "spearman")
  
  corr_results <- rbind(corr_results, data.frame(
    gene = gene,
    rho = round(ct$estimate, 3),
    pval = ct$p.value,
    n_probes = length(probes)
  ))
}

corr_results$padj <- p.adjust(corr_results$pval, method = "BH")
corr_results <- corr_results[order(corr_results$padj), ]

write.csv(corr_results, "03_results/methylation_expression_correlation.csv", row.names = FALSE)
cat("Correlation results:\n")
print(corr_results)