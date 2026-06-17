# 22_methylation_figures.R
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)

# Load methylation results
results_meth <- readRDS("03_results/dmeth_tcga_tumor_vs_normal.rds")

# ── Volcano Plot ──────────────────────────────────────────
results_meth$sig <- ifelse(results_meth$adj.P.Val < 0.05 & 
                             results_meth$logFC > 0.1, "Hypermethylated",
                           ifelse(results_meth$adj.P.Val < 0.05 & 
                                    results_meth$logFC < -0.1, "Hypomethylated", 
                                  "Not significant"))

# Load nexus probe annotation
meth_se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
probe_annot <- as.data.frame(rowData(meth_se))
probe_annot$probe_id <- rownames(probe_annot)

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

nexus_probes <- probe_annot[probe_annot$gene %in% nexus_genes, 
                            c("probe_id", "gene")]

# Add gene labels to results
results_meth$probe_id <- rownames(results_meth)
results_meth$gene_label <- nexus_probes$gene[match(results_meth$probe_id, 
                                                   nexus_probes$probe_id)]

# Volcano
p_volc <- ggplot(results_meth, aes(x = logFC, y = -log10(adj.P.Val),
                                   color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_text_repel(data = results_meth[!is.na(results_meth$gene_label) & 
                                        results_meth$adj.P.Val < 0.05, ],
                  aes(label = gene_label), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Hypermethylated" = "#B2182B",
                                "Hypomethylated" = "#2166AC",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "Differential Methylation — TCGA-LIHC (Tumor vs Normal)",
       x = "Delta Beta (logFC)", y = "-log10(adj.P.Val)",
       color = "Status") +
  theme_bw()

png("04_figures/volcano_methylation.png", width = 10, height = 7, 
    units = "in", res = 300)
print(p_volc)
dev.off()
cat("Methylation volcano saved.\n")


# ── Nexus Gene Methylation Heatmap ───────────────────────
library(ComplexHeatmap)
library(circlize)

# Load beta matrix and filter to nexus probes
beta_mat <- readRDS("data/raw/tcga_lihc_meth_filtered.rds")

# Keep only nexus probes
nexus_probe_ids <- nexus_probes$probe_id[nexus_probes$probe_id %in% rownames(beta_mat)]
beta_nexus <- beta_mat[nexus_probe_ids, ]

# Add gene labels as row names
rownames(beta_nexus) <- nexus_probes$gene[match(rownames(beta_nexus), 
                                                nexus_probes$probe_id)]

# Sample annotation
keep <- meth_se$shortLetterCode %in% c("TP", "NT")
sample_type <- meth_se$shortLetterCode[keep]

ha <- HeatmapAnnotation(
  Type = sample_type,
  col = list(Type = c("NT" = "#2E7C4A", "TP" = "#D9534F")),
  show_legend = TRUE
)

# Plot
png("04_figures/heatmap_methylation_nexus.png", 
    width = 14, height = 10, units = "in", res = 300)
Heatmap(beta_nexus,
        name = "Beta value",
        top_annotation = ha,
        show_column_names = FALSE,
        col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "white", "#B2182B")),
        column_title = "Nexus Gene Methylation — TCGA-LIHC",
        cluster_rows = TRUE,
        cluster_columns = TRUE)
dev.off()
cat("Methylation heatmap saved.\n")


png("04_figures/heatmap_methylation_nexus.png", 
    width = 14, height = 10, units = "in", res = 300)
Heatmap(beta_nexus,
        name = "Beta value",
        top_annotation = ha,
        show_column_names = FALSE,
        show_row_names = FALSE,
        col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "white", "#B2182B")),
        column_title = "Nexus Gene Methylation — TCGA-LIHC",
        cluster_rows = TRUE,
        cluster_columns = TRUE)
dev.off()