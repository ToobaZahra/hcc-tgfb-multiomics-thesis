# 04_normalize_pca.R
library(DESeq2)

dds <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")

# VST normalization (needed for PCA + heatmap)
vsd <- vst(dds, blind = FALSE)

# Save vsd object
saveRDS(vsd, "data/raw/tcga_lihc_vsd.rds")
cat("VST normalization done. Saved to data/raw/tcga_lihc_vsd.rds\n")

# PCA plot
library(ggplot2)
pcaData <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: TCGA-LIHC Tumor vs Normal") +
  theme_bw()

ggsave("04_figures/pca_tumor_vs_normal.png", p, width = 8, height = 6, dpi = 300)
cat("PCA plot saved.\n")