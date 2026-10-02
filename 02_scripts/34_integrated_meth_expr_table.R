library(dplyr); library(ggplot2)
library(org.Hs.eg.db)
select <- dplyr::select

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")

# ── 1. Methylation + correlation (from script 33) ────────
summ <- read.csv("03_results/methylation/meth_promoter_gene_summary.csv")
cor_df <- read.csv("03_results/methylation/meth_expr_corr_promoter_tumor.csv")

# ── 2. DGE tumor vs normal (from script 06) ──────────────
dge <- readRDS("03_results/dge/dge_tcga_tumor_vs_normal.rds")
dge <- as.data.frame(dge)
cat("DGE columns:", paste(colnames(dge), collapse = ", "), "\n")

sym_col <- intersect(c("symbol","gene_name","gene","SYMBOL","external_gene_name"), colnames(dge))
if (length(sym_col) > 0) {
  dge$gene <- dge[[sym_col[1]]]
} else {
  ens <- sub("\\..*$", "", rownames(dge))
  dge$gene <- suppressMessages(mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first"))
}
nexus_genes <- summ$gene
dge_n <- dge %>% filter(gene %in% nexus_genes) %>%
  group_by(gene) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(gene, log2FC = log2FoldChange, dge_padj = padj)
cat("Nexus genes found in DGE:", nrow(dge_n), "of", length(nexus_genes), "\n")

# ── 3. Join + classify ───────────────────────────────────
tab <- summ %>%
  select(gene, n_probes, beta_normal, beta_tumor, delta, pct_hyper, pct_hypo) %>%
  left_join(cor_df %>% select(gene, rho, cor_padj = padj), by = "gene") %>%
  left_join(dge_n, by = "gene") %>%
  mutate(
    meth_status = case_when(delta >= 0.05 ~ "Hyper", delta <= -0.05 ~ "Hypo", TRUE ~ "Stable"),
    expr_status = case_when(dge_padj < 0.05 & log2FC >= 1 ~ "Up",
                            dge_padj < 0.05 & log2FC <= -1 ~ "Down", TRUE ~ "NS"),
    link = case_when(cor_padj < 0.05 & rho < 0 ~ "Negative",
                     cor_padj < 0.05 & rho > 0 ~ "Positive", TRUE ~ "None"),
    interpretation = case_when(
      meth_status == "Hyper" & link == "Negative" & expr_status == "Down" ~ "Silenced by hypermethylation",
      meth_status == "Hypo"  & link == "Negative" & expr_status == "Up"   ~ "Activated by hypomethylation",
      meth_status != "Stable" & link == "Negative" ~ "Link, but tumor expression not concordant",
      meth_status != "Stable" & link == "Positive" ~ "Unexpected direction",
      meth_status != "Stable" & link == "None" ~ "Methylation shift, no expression link",
      meth_status == "Stable" & link == "Negative" ~ "Expression link, no tumor methylation shift",
      TRUE ~ "No clear effect")
  ) %>% arrange(desc(delta))

write.csv(tab, "03_results/methylation/nexus_integrated_meth_expr_table.csv", row.names = FALSE)
print(tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        select(gene, delta, meth_status, rho, link, log2FC, expr_status, interpretation))
print(table(tab$interpretation))

# ── 4. Figure: methylation shift vs expression change ────
cols <- c("Silenced by hypermethylation" = "#D7301F",
          "Activated by hypomethylation" = "#2C7BB6",
          "Unexpected direction" = "#FF7F00",
          "Methylation shift, no expression link" = "#6A3D9A",
          "Expression link, no tumor methylation shift" = "#33A02C",
          "No clear effect" = "grey60",
          "Link, but tumor expression not concordant" = "#E7298A")
p <- ggplot(tab, aes(delta, log2FC, color = interpretation, size = abs(rho))) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = cols, name = NULL) +
  scale_size_continuous(range = c(2, 6), name = "|rho|") +
  labs(title = "Promoter methylation change vs expression change in HCC",
       subtitle = "x = tumor - normal promoter beta; y = log2 fold change (DESeq2); colour = methylation-expression link",
       x = "Delta beta (tumor - normal)", y = "log2 fold change (tumor vs normal)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.6, show.legend = FALSE)
} else {
  p <- p + geom_text(aes(label = gene), size = 3, vjust = -1, show.legend = FALSE)
}
ggsave("04_figures/fig_meth_vs_dge_quadrant.png", p, width = 11, height = 7.5, dpi = 300, bg = "white")
print(p)
