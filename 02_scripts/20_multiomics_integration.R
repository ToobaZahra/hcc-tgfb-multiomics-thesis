# 20_multiomics_integration.R
library(dplyr)

# Load all results
dge <- read.csv("03_results/dge_nexus_genes.csv")
meth <- read.csv("03_results/dmeth_nexus_summary.csv")
cox <- read.csv("03_results/cox_nexus_genes.csv")
replication <- read.csv("03_results/nexus_replication_summary.csv")

# Clean DGE
dge_clean <- dge %>%
  select(symbol, logFC_TCGA = log2FoldChange, padj_TCGA = padj) %>%
  mutate(DE_sig = padj_TCGA < 0.05,
         DE_dir = ifelse(logFC_TCGA > 0, "UP", "DOWN"))



# Remove NA rows
integrated_clean <- integrated[!is.na(integrated$symbol), ]

write.csv(integrated_clean, "03_results/multiomics_integrated.csv", row.names = FALSE)

cat("Final multi-omics hits (score >= 3):\n")
print(integrated_clean[integrated_clean$omics_score >= 3,
                       c("symbol", "DE_dir", "meth_dir", "cox_dir",
                         "replicated", "omics_score")])

# Clean methylation
meth_clean <- meth %>%
  select(symbol = gene, n_sig_probes = n_sig,
         n_hyper, n_hypo, mean_logFC_meth = mean_logFC) %>%
  mutate(meth_sig = n_sig_probes > 0,
         meth_dir = ifelse(mean_logFC_meth > 0, "HYPER", "HYPO"))

# Clean Cox
cox_clean <- cox %>%
  select(symbol = gene, HR, pval_cox = pval) %>%
  mutate(cox_sig = pval_cox < 0.05,
         cox_dir = ifelse(HR > 1, "RISK", "PROTECTIVE"))

# Clean replication
rep_clean <- replication %>%
  select(symbol, TCGA_dir, GSE14520_dir, GSE36376_dir,
         GSE76427_dir, replicated = consistent)

# Merge all
integrated <- dge_clean %>%
  left_join(meth_clean, by = "symbol") %>%
  left_join(cox_clean, by = "symbol") %>%
  left_join(rep_clean, by = "symbol")

# Score each gene (how many omics layers significant)
integrated <- integrated %>%
  mutate(omics_score = as.integer(DE_sig) +
           as.integer(meth_sig) +
           as.integer(cox_sig) +
           as.integer(replicated == TRUE))

integrated <- integrated %>% arrange(desc(omics_score), padj_TCGA)

integrated_clean <- integrated[!is.na(integrated$symbol) & 
                                 !is.na(integrated$omics_score), ]

write.csv(integrated_clean, "03_results/multiomics_integrated.csv", row.names = FALSE)

cat("Final multi-omics hits (score >= 3):\n")
print(integrated_clean[integrated_clean$omics_score >= 3,
                       c("symbol", "DE_dir", "meth_dir", "cox_dir",
                         "replicated", "omics_score")])


library(ggplot2)

# Build plot data
plot_df <- integrated_clean[integrated_clean$omics_score >= 3, ]

# Reshape for heatmap
library(tidyr)
heat_df <- data.frame(
  symbol = rep(plot_df$symbol, 4),
  layer  = c(rep("Expression", 7), rep("Methylation", 7),
             rep("Survival HR", 7), rep("Replicated", 7)),
  value  = c(ifelse(plot_df$DE_dir == "UP", 1, -1),
             ifelse(plot_df$meth_dir == "HYPER", 1, -1),
             ifelse(plot_df$cox_dir == "RISK", 1, -1),
             ifelse(plot_df$replicated == TRUE, 1, 0))
)

p <- ggplot(heat_df, aes(x = layer, y = symbol, fill = factor(value))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("-1" = "#2166AC", "0" = "grey80", "1" = "#B2182B"),
                    labels = c("Down/Hypo/Protective/No", "Not replicated", 
                               "Up/Hyper/Risk/Yes")) +
  labs(title = "Multi-omics Integration — TGF-β Nexus Genes",
       x = "", y = "", fill = "Direction") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        plot.title = element_text(size = 13))

png("04_figures/multiomics_heatmap.png", width = 8, height = 5, 
    units = "in", res = 300)
print(p)
dev.off()
cat("Heatmap saved.\n")