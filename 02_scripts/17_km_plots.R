# 17_km_plots.R
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)

# Load matched data
clin_tumor <- readRDS("03_results/survival_clinical_matched.rds")
vsd_tumor <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds")

# Load DEG results for gene symbols
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

# Nexus genes to plot
nexus_genes <- c("TGFB1","TGFBR2","SMAD3","SMAD4","MYC","KLF4",
                 "DUSP1","HDAC11","SNAI1","IL6")

# Get expression matrix
expr_mat <- assay(vsd_tumor)

# Match gene symbols to Ensembl IDs
gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

# Create output folder
dir.create("04_figures/km_plots", showWarnings = FALSE, recursive = TRUE)

# Loop through nexus genes and plot KM
for (i in 1:nrow(gene_map)) {
  sym <- gene_map$symbol[i]
  ens <- gene_map$ensembl_id[i]
  
  # Find matching row in expr_mat
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) == 0) next
  
  expr <- expr_mat[row_match[1], ]
  median_val <- median(expr)
  group <- ifelse(expr >= median_val, "High", "Low")
  
  surv_df <- data.frame(
    os_time  = clin_tumor$os_time,
    os_event = clin_tumor$os_event,
    group    = factor(group, levels = c("Low", "High"))
  )
  surv_df <- surv_df[!is.na(surv_df$os_time), ]
  
  fit <- survfit(Surv(os_time, os_event) ~ group, data = surv_df)
  
  p <- ggsurvplot(fit,
                  data = surv_df,
                  pval = TRUE,
                  risk.table = TRUE,
                  title = paste0(sym, " — OS (TCGA-LIHC)"),
                  legend.labs = c("Low", "High"),
                  palette = c("#2E7C4A", "#D9534F"))
  
  png(paste0("04_figures/km_plots/km_", sym, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
  cat("Saved KM plot for", sym, "\n")
}

cat("All KM plots done.\n")