## 35: probe-level promoter methylation vs expression
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/methylation"; fig_dir <- "04_figures"
genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# 1. Promoter probes
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit  <- vapply(g_list, function(g) any(g %in% genes), logical(1))
keep <- hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)
cand <- ann[keep, ]; cand$gene_list <- g_list[keep]

ent <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL"))
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
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]

# 2. Distance to TSS
tss_tab <- data.frame(gene   = prom$symbol,
                      strand = as.character(strand(prom)),
                      tss    = ifelse(as.character(strand(prom)) == "+",
                                      start(prom) + 1500, end(prom) - 1500)) %>% distinct()
probe_gene$pos <- cand$beg_A[match(probe_gene$probe, cand$probe)]
probe_gene$dist_TSS <- mapply(function(g, p) {
  t <- tss_tab[tss_tab$gene == g, ]
  d <- ifelse(t$strand == "+", p - t$tss, t$tss - p)
  d[which.min(abs(d))]
}, probe_gene$gene, probe_gene$pos)
cat("Promoter probes per gene:\n"); print(table(probe_gene$gene))

# 3. Tumor beta per probe
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
beta_t <- assay(se[unique(probe_gene$probe), tum])
colnames(beta_t) <- substr(colnames(beta_t), 1, 12)

# 4. Tumor expression
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
idx <- match(sym2ens, sub("\\..*$", "", rownames(expr)))
ok  <- !is.na(idx)
expr <- expr[idx[ok], , drop = FALSE]; rownames(expr) <- names(sym2ens)[ok]
et <- which(substr(colnames(expr), 14, 15) == "01")
et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
expr_t <- expr[, et, drop = FALSE]; colnames(expr_t) <- substr(colnames(expr_t), 1, 12)
common <- intersect(colnames(beta_t), colnames(expr_t))
cat("Patients with methylation + expression:", length(common), "\n")

# 5. Spearman per probe (skips probes with too little data)
res <- do.call(rbind, lapply(seq_len(nrow(probe_gene)), function(i) {
  g <- probe_gene$gene[i]; p <- probe_gene$probe[i]
  if (!g %in% rownames(expr_t)) return(NULL)
  x <- beta_t[p, common]; y <- expr_t[g, common]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 30 || sd(x[ok]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.frame(gene = g, probe = p, dist_TSS = probe_gene$dist_TSS[i],
             mean_beta = mean(x[ok]), n = sum(ok),
             rho = unname(ct$estimate), p = ct$p.value)
}))
cat("Probes tested:", nrow(res), "of", nrow(probe_gene), "\n")
res$padj <- p.adjust(res$p, "BH")
res <- res %>% arrange(gene, rho)
write.csv(res, file.path(out_dir, "meth_expr_corr_probe_level.csv"), row.names = FALSE)

cat("\nBest (most negative rho) probe per gene:\n")
print(res %>% group_by(gene) %>% slice_min(rho, n = 1) %>%
        select(gene, probe, dist_TSS, mean_beta, rho, padj))

# 6. Plot
p <- ggplot(res, aes(dist_TSS, rho, color = padj < 0.05)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) +
  facet_wrap(~ gene, ncol = 1) +
  scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D7301F"), name = "FDR < 0.05") +
  labs(title = "Probe-level promoter methylation vs expression",
       subtitle = "Each dot = one CpG probe; dashed line = TSS; negative rho = more methylation, less expression",
       x = "Distance to TSS (bp)", y = "Spearman rho") +
  theme_minimal(base_size = 12)
ggsave(file.path(fig_dir, "fig_meth_expr_probe_level.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")

r <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
sel <- r %>% filter((rho < -0.3 | rho > 0.3) & padj < 0.05) %>% arrange(rho)
write.csv(sel, "03_results/methylation/probes_rho_beyond_0.3.csv", row.names = FALSE)
print(as.data.frame(sel[, c("gene","probe","dist_TSS","mean_beta","rho","padj")]))