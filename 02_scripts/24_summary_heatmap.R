# 24_summary_heatmap.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

integrated <- read.csv("03_results/multiomics_integrated.csv")

plot_df <- integrated[!is.na(integrated$omics_score) & 
                        integrated$omics_score >= 2, ]

heat_data <- plot_df %>%
  select(symbol, logFC_TCGA, mean_logFC_meth, HR, omics_score) %>%
  column_to_rownames("symbol")

heat_scaled <- as.data.frame(scale(heat_data))

heat_long <- heat_scaled %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "metric", values_to = "value")

p <- ggplot(heat_long, aes(x = metric, y = gene, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Z-score") +
  labs(title = "Multi-Omics Summary Heatmap — TGF-β Nexus Genes",
       x = "", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("04_figures/summary_heatmap.png", width = 8, height = 6,
    units = "in", res = 300)
print(p)
dev.off()
cat("Summary heatmap saved.\n")