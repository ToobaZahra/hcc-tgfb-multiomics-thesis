# 06_dge_tcga.R
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

dds <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")

# Run DESeq (slow: 10–30 minutes)
dds <- DESeq(dds)

# Extract results: tumor vs normal
res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
summary(res)

# Annotate
res_df <- as.data.frame(res)
res_df$ensembl_id <- gsub("\\..*$", "", rownames(res_df))
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = res_df$ensembl_id,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

dir.create("03_results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res_df, "03_results/dge_tcga_tumor_vs_normal.rds")
write.csv(res_df, "03_results/dge_tcga_tumor_vs_normal.csv", row.names = FALSE)

# Save full dds for later
saveRDS(dds, "data/raw/tcga_lihc_dds_DEseqed.rds")


nexus_genes <- c("TGFB1","TGFBR1","TGFBR2",
                 "SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10",
                 "DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1",
                 "TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

nexus <- res_df[res_df$symbol %in% nexus_genes, ]
nexus <- nexus[order(nexus$padj), ]
write.csv(nexus, "03_results/dge_nexus_genes.csv", row.names = FALSE)
print(nexus[, c("symbol", "log2FoldChange", "padj")])