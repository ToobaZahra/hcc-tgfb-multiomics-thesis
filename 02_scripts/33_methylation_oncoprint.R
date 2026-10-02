library(SummarizedExperiment)
library(DESeq2)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/methylation"; fig_dir <- "04_figures"
dir.create(fig_dir, showWarnings = FALSE)
delta_thr <- 0.10    # |tumor beta - normal mean beta| for hyper/hypo call

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Probes -> nexus genes (exact match, multi-gene safe) ─
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit <- vapply(g_list, function(g) any(g %in% nexus_genes), logical(1))
cand <- ann[hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A), ]
cand$gene_list <- g_list[hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)]
cat("Candidate probes touching nexus genes:", nrow(cand), "\n")

# ── 2. Promoter windows from TSS (-1500/+500) ────────────
ent <- suppressMessages(mapIds(org.Hs.eg.db, nexus_genes, "ENTREZID", "SYMBOL"))
cat("Genes without Entrez ID:", paste(names(ent)[is.na(ent)], collapse = ", "), "\n")
ent <- ent[!is.na(ent)]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, "gene")[ent]
prom <- unlist(promoters(tx_by_gene, upstream = 1500, downstream = 500))
prom$symbol <- names(ent)[match(names(prom), ent)]
prom <- keepStandardChromosomes(prom, pruning.mode = "coarse")

probe_gr <- GRanges(cand$chrm_A, IRanges(cand$beg_A, cand$beg_A + 1))
ov <- suppressWarnings(findOverlaps(probe_gr, prom))
probe_gene <- unique(data.frame(probe = cand$probe[queryHits(ov)],
                                gene  = prom$symbol[subjectHits(ov)]))
# keep pair only if the probe is annotated to that gene
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]
cat("\nPromoter probes per gene:\n"); print(table(factor(probe_gene$gene, levels = nexus_genes)))

# ── 3. Gene-level promoter beta (mean of probes) ─────────
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
beta <- assay(se[unique(probe_gene$probe), c(nor, tum)])

gene_beta <- do.call(rbind, lapply(split(probe_gene$probe, probe_gene$gene),
                                   function(p) colMeans(beta[p, , drop = FALSE], na.rm = TRUE)))
n_nor <- length(nor)
b_nor <- gene_beta[, seq_len(n_nor), drop = FALSE]
b_tum <- gene_beta[, -seq_len(n_nor), drop = FALSE]
colnames(b_tum) <- substr(colnames(b_tum), 1, 12)

# ── 4. Tumor vs normal per gene ──────────────────────────
summ <- data.frame(gene = rownames(gene_beta),
                   n_probes = as.integer(table(probe_gene$gene)[rownames(gene_beta)]),
                   beta_normal = rowMeans(b_nor, na.rm = TRUE),
                   beta_tumor  = rowMeans(b_tum, na.rm = TRUE))
summ$delta <- summ$beta_tumor - summ$beta_normal
summ$p <- sapply(seq_len(nrow(b_tum)), function(i)
  wilcox.test(b_tum[i, ], b_nor[i, ])$p.value)
summ$padj <- p.adjust(summ$p, "BH")

delta_mat <- b_tum - rowMeans(b_nor, na.rm = TRUE)          # gene x patient
calls <- ifelse(delta_mat > delta_thr, "Hyper",
                ifelse(delta_mat < -delta_thr, "Hypo", "Neutral"))
summ$pct_hyper <- 100 * rowMeans(calls == "Hyper", na.rm = TRUE)
summ$pct_hypo  <- 100 * rowMeans(calls == "Hypo",  na.rm = TRUE)
summ <- summ %>% arrange(desc(delta))
write.csv(summ, file.path(out_dir, "meth_promoter_gene_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(t(calls)) %>% tibble::rownames_to_column("patient"),
          file.path(out_dir, "meth_patient_calls.csv"), row.names = FALSE)
print(summ %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ── 5. Oncoprints ────────────────────────────────────────
plot_onco <- function(genes, file, title) {
  cl <- calls[genes, , drop = FALSE]
  pat_ord <- names(sort(colSums(cl == "Hyper", na.rm = TRUE) -
                          colSums(cl == "Hypo",  na.rm = TRUE), decreasing = TRUE))
  lab <- setNames(sprintf("%s  (%.0f%% hyper | %.0f%% hypo)", genes,
                          summ$pct_hyper[match(genes, summ$gene)],
                          summ$pct_hypo[match(genes, summ$gene)]), genes)
  df <- as.data.frame(as.table(cl), stringsAsFactors = FALSE)
  names(df) <- c("gene", "patient", "call")
  df$gene <- factor(lab[df$gene], levels = rev(lab[genes]))
  df$patient <- factor(df$patient, levels = pat_ord)
  p <- ggplot(df, aes(patient, gene, fill = call)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = c(Hyper = "#D7301F", Neutral = "grey92", Hypo = "#2C7BB6"),
                      na.value = "white", name = NULL) +
    labs(title = title,
         subtitle = sprintf("TCGA-LIHC tumors (n = %d); call = promoter beta vs normal mean, |delta| > %.2f",
                            length(pat_ord), delta_thr),
         x = paste0(length(pat_ord), " patients"), y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom")
  ggsave(file, p, width = 11, height = max(3.5, 0.28 * length(genes) + 2), dpi = 300, bg = "white")
  print(p)
}
plot_onco(summ$gene, file.path(fig_dir, "fig_meth_oncoprint_all.png"),
          "Promoter methylation oncoprint: nexus genes (hyper to hypo)")
top <- c(head(summ$gene, 5), tail(summ$gene, 5))
plot_onco(top, file.path(fig_dir, "fig_meth_oncoprint_top.png"),
          "Top 5 hypermethylated and top 5 hypomethylated nexus genes")

# ── 6. Methylation vs expression (tumor only, matched patients) ─
library(DESeq2)
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds")
expr <- assay(vsd)

# map nexus symbols -> Ensembl IDs, then match to vsd rows (version suffix removed)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, nexus_genes, "ENSEMBL", "SYMBOL",
                                   multiVals = "first"))
cat("Genes with no Ensembl ID:", paste(names(sym2ens)[is.na(sym2ens)], collapse = ", "), "\n")
row_ens <- sub("\\..*$", "", rownames(expr))
idx <- match(sym2ens, row_ens)
cat("Genes not found in vsd:", paste(names(sym2ens)[is.na(idx)], collapse = ", "), "\n")
ok <- !is.na(idx)
expr <- expr[idx[ok], , drop = FALSE]
rownames(expr) <- names(sym2ens)[ok]

et <- which(substr(colnames(expr), 14, 15) == "01")
et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
expr_t <- expr[, et, drop = FALSE]; colnames(expr_t) <- substr(colnames(expr_t), 1, 12)

common <- intersect(colnames(b_tum), colnames(expr_t))
cat("\nPatients with methylation + expression:", length(common), "\n")
cg <- intersect(rownames(b_tum), rownames(expr_t))
cor_df <- do.call(rbind, lapply(cg, function(g) {
  ct <- suppressWarnings(cor.test(b_tum[g, common], expr_t[g, common], method = "spearman"))
  data.frame(gene = g, rho = unname(ct$estimate), p = ct$p.value)
}))
cor_df$padj <- p.adjust(cor_df$p, "BH")
cor_df <- cor_df %>% left_join(summ[, c("gene", "n_probes", "delta")], by = "gene") %>%
  arrange(rho)
write.csv(cor_df, file.path(out_dir, "meth_expr_corr_promoter_tumor.csv"), row.names = FALSE)
print(cor_df %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

p_cor <- ggplot(cor_df, aes(reorder(gene, rho), rho, fill = padj < 0.05)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#2C7BB6", `FALSE` = "grey70"), name = "FDR < 0.05") +
  labs(title = "Promoter methylation vs expression (Spearman)",
       subtitle = "Tumor samples only; negative rho = higher methylation, lower expression",
       x = NULL, y = "Spearman rho") +
  theme_minimal(base_size = 12) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "fig_meth_expr_rho.png"), p_cor, width = 7, height = 6, dpi = 300, bg = "white")

sel <- head(cor_df$gene[order(-abs(cor_df$rho))], 6)
sc_df <- do.call(rbind, lapply(sel, function(g)
  data.frame(gene = g, beta = b_tum[g, common], expr = expr_t[g, common])))
sc_df$gene <- factor(sc_df$gene, levels = sel)
p_sc <- ggplot(sc_df, aes(beta, expr)) +
  geom_point(alpha = 0.5, size = 1.2) + geom_smooth(method = "lm", se = FALSE, color = "#D7301F") +
  facet_wrap(~ gene, scales = "free") +
  labs(title = "Promoter beta vs expression: 6 strongest genes",
       x = "Promoter mean beta", y = "VST expression") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_meth_expr_scatter.png"), p_sc, width = 9, height = 6, dpi = 300, bg = "white")
print(p_cor); print(p_sc)