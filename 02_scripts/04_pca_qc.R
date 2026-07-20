# 04_pca_qc.R
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load filtered DDS
dds <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")

# VST normalization
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "data/raw/tcga_lihc_vsd.rds")
cat("VST done. Dimensions:", dim(vsd), "\n")

# ── PCA ──────────────────────────────────────────────────
pcaData <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 2, alpha = 0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("Primary Tumor" = "#D9534F",
                                "Solid Tissue Normal" = "#2E7C4A")) +
  ggtitle("PCA — TCGA-LIHC (VST normalized)") +
  theme_bw()

png("04_figures/pca_tcga_lihc.png", width = 8, height = 6, units = "in", res = 300)
print(p_pca)
dev.off()
cat("PCA plot saved.\n")

# ── Sample distance heatmap ───────────────────────────────
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Annotation
annot <- data.frame(
  Type = vsd$sample_type,
  row.names = colnames(vsd)
)

png("04_figures/sample_dist_heatmap.png", width = 12, height = 10,
    units = "in", res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_row = annot,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Distance Heatmap — TCGA-LIHC")
dev.off()
cat("Sample distance heatmap saved.\n")