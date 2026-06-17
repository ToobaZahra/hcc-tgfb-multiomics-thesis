# 19_pathway_enrichment.R
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# Load DEG results
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

# Filter significant DEGs
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                abs(res_df$log2FoldChange) > 1, ]
cat("Significant DEGs:", nrow(sig), "\n")

# Convert Ensembl to Entrez IDs
sig$entrez <- mapIds(org.Hs.eg.db,
                     keys = sig$ensembl_id,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

sig <- sig[!is.na(sig$entrez), ]
cat("With Entrez IDs:", nrow(sig), "\n")

# GO Biological Process enrichment
ego <- enrichGO(gene          = sig$entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

cat("GO-BP terms enriched:", nrow(ego@result[ego@result$p.adjust < 0.05, ]), "\n")

# KEGG pathway enrichment
ekegg <- enrichKEGG(gene         = sig$entrez,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

cat("KEGG pathways enriched:", nrow(ekegg@result[ekegg@result$p.adjust < 0.05, ]), "\n")

# Save results
write.csv(ego@result, "03_results/enrichment_GO_BP.csv", row.names = FALSE)
write.csv(ekegg@result, "03_results/enrichment_KEGG.csv", row.names = FALSE)

# Dotplot GO
p1 <- dotplot(ego, showCategory = 20, title = "GO Biological Process") +
  theme(axis.text.y = element_text(size = 8))

png("04_figures/dotplot_GO_BP.png", width = 10, height = 8, units = "in", res = 300)
print(p1)
dev.off()

# Dotplot KEGG
p2 <- dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment") +
  theme(axis.text.y = element_text(size = 8))

png("04_figures/dotplot_KEGG.png", width = 10, height = 8, units = "in", res = 300)
print(p2)
dev.off()

cat("Plots saved.\n")



# GSEA on TGF-β related terms
tgfb_terms <- ego@result[grep("TGF|transform|smad|fibrosis|epithelial", 
                              ego@result$Description, 
                              ignore.case = TRUE), ]

cat("\nTGF-β related GO terms:\n")
print(tgfb_terms[, c("Description", "p.adjust", "Count")])

write.csv(tgfb_terms, "03_results/enrichment_TGFB_terms.csv", row.names = FALSE)


# GSEA with ranked gene list
library(clusterProfiler)

# Create ranked list — all genes sorted by log2FC
res_df_clean <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$ensembl_id), ]

# Convert to Entrez
res_df_clean$entrez <- mapIds(org.Hs.eg.db,
                              keys = res_df_clean$ensembl_id,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

res_df_clean <- res_df_clean[!is.na(res_df_clean$entrez), ]

# Ranked vector
gene_list <- res_df_clean$log2FoldChange
names(gene_list) <- res_df_clean$entrez
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

cat("Ranked gene list length:", length(gene_list), "\n")

# GSEA GO-BP
gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

cat("GSEA GO-BP terms:", nrow(gsea_go@result[gsea_go@result$p.adjust < 0.05, ]), "\n")

# Save
write.csv(gsea_go@result, "03_results/gsea_GO_BP.csv", row.names = FALSE)

# Plot top 20
p3 <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(axis.text.y = element_text(size = 7))

png("04_figures/gsea_GO_BP.png", width = 14, height = 8, units = "in", res = 300)
print(p3)
dev.off()
cat("GSEA plot saved.\n")