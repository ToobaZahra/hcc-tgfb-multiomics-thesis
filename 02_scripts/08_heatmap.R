# 08_heatmap.R
library(ComplexHeatmap)
library(circlize)
library(SummarizedExperiment)

vsd <- readRDS("data/raw/tcga_lihc_vsd.rds")
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

# Top 50 by padj (excluding NA)
res_df <- res_df[!is.na(res_df$padj), ]
top50 <- head(res_df[order(res_df$padj), ], 50)

mat <- assay(vsd)[rownames(top50), ]
mat <- t(scale(t(mat)))   # z-score by gene

ha <- HeatmapAnnotation(
  Type = vsd$sample_type,
  col = list(Type = c("Solid Tissue Normal" = "#2E7C4A",
                      "Primary Tumor" = "#D9534F")),
  show_legend = TRUE
)

png("04_figures/heatmap_top50.png", width = 12, height = 10, units = "in", res = 300)
Heatmap(mat,
        name = "z-score",
        top_annotation = ha,
        show_column_names = FALSE,
        row_labels = top50$symbol,
        col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
        column_title = "Top 50 DEGs (TCGA-LIHC)")
dev.off()
cat("Heatmap saved.\n")