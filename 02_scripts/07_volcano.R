# 07_volcano.R
library(EnhancedVolcano)
library(ggplot2)

res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1")

p <- EnhancedVolcano(
  res_df,
  lab = res_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  selectLab = nexus_genes,
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TCGA-LIHC: Tumor vs Normal",
  subtitle = "TGF-β nexus genes labelled",
  caption = paste0("Total genes: ", nrow(res_df)),
  pointSize = 1.5,
  labSize = 3.5,
  drawConnectors = TRUE,
  widthConnectors = 0.4
)

ggsave("04_figures/volcano_tcga.png", p, width = 9, height = 7, dpi = 300)
ggsave("04_figures/volcano_tcga.pdf", p, width = 9, height = 7)
cat("Volcano plot saved.\n")