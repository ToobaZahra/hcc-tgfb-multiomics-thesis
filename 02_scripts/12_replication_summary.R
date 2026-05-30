# 12_replication_summary.R
# Summarize nexus gene replication across TCGA + 3 GEO cohorts

library(dplyr)

# Load all nexus results
tcga <- read.csv("03_results/dge_nexus_genes.csv")
gse14520 <- read.csv("03_results/limma_GSE14520_nexus_genes.csv")
gse36376 <- read.csv("03_results/limma_GSE36376_nexus_genes.csv")
gse76427 <- read.csv("03_results/limma_GSE76427_nexus_genes.csv")

# Helper: get best probe per gene (lowest adj.P.Val)
best_probe <- function(df, fc_col, p_col) {
  df %>%
    group_by(symbol) %>%
    slice_min(order_by = .data[[p_col]], n = 1) %>%
    ungroup() %>%
    select(symbol, logFC = all_of(fc_col), adj.P.Val = all_of(p_col))
}

tcga_clean    <- tcga %>% select(symbol, logFC = log2FoldChange, adj.P.Val = padj)
gse14520_clean <- best_probe(gse14520, "logFC", "adj.P.Val")
gse36376_clean <- best_probe(gse36376, "logFC", "adj.P.Val")
gse76427_clean <- best_probe(gse76427, "logFC", "adj.P.Val")

# Merge all
summary_df <- tcga_clean %>%
  rename(TCGA_logFC = logFC, TCGA_padj = adj.P.Val) %>%
  left_join(gse14520_clean %>% rename(GSE14520_logFC = logFC, GSE14520_padj = adj.P.Val), by = "symbol") %>%
  left_join(gse36376_clean %>% rename(GSE36376_logFC = logFC, GSE36376_padj = adj.P.Val), by = "symbol") %>%
  left_join(gse76427_clean %>% rename(GSE76427_logFC = logFC, GSE76427_padj = adj.P.Val), by = "symbol")

# Add direction consistency
summary_df <- summary_df %>%
  mutate(
    TCGA_dir    = ifelse(TCGA_logFC > 0, "UP", "DOWN"),
    GSE14520_dir = ifelse(GSE14520_logFC > 0, "UP", "DOWN"),
    GSE36376_dir = ifelse(GSE36376_logFC > 0, "UP", "DOWN"),
    GSE76427_dir = ifelse(GSE76427_logFC > 0, "UP", "DOWN"),
    consistent = (TCGA_dir == GSE14520_dir) & 
      (TCGA_dir == GSE36376_dir) & 
      (TCGA_dir == GSE76427_dir)
  )

write.csv(summary_df, "03_results/nexus_replication_summary.csv", row.names = FALSE)

cat("Replication summary:\n")
print(summary_df[, c("symbol", "TCGA_dir", "GSE14520_dir", 
                     "GSE36376_dir", "GSE76427_dir", "consistent")])
cat("\nConsistent across all 4 cohorts:", sum(summary_df$consistent, na.rm = TRUE), "genes\n")