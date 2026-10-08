## 41: association of methylation with grade, G1-G4 as an ordered variable
## (all promoter probes of the 26 genes + gene-level mean beta), tumors only
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/methylation"; fig_dir <- "04_figures"

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Promoter probes (same logic as scripts 33/35/39) ──
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit  <- vapply(g_list, function(g) any(g %in% genes), logical(1))
keep <- hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)
cand <- ann[keep, ]; cand$gene_list <- g_list[keep]

ent <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")); ent <- ent[!is.na(ent)]
tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")[ent]
prom <- unlist(promoters(tx_by_gene, upstream = 1500, downstream = 500))
prom$symbol <- names(ent)[match(names(prom), ent)]
prom <- keepStandardChromosomes(prom, pruning.mode = "coarse")
probe_gr <- GRanges(cand$chrm_A, IRanges(cand$beg_A, cand$beg_A + 1))
ov <- suppressWarnings(findOverlaps(probe_gr, prom))
probe_gene <- unique(data.frame(probe = cand$probe[queryHits(ov)], gene = prom$symbol[subjectHits(ov)]))
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]
cat("Promoter probes:", length(unique(probe_gene$probe)), "\n")

# ── 2. Tumors + grade as an ordered number (G1=1 ... G4=4) ──
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cd  <- as.data.frame(colData(se))[tum, ]
cd$grade_num <- match(cd$tumor_grade, c("G1", "G2", "G3", "G4"))      # GX / missing -> NA
cd$grade_grp <- ifelse(cd$grade_num %in% 1:2, "Low", ifelse(cd$grade_num %in% 3:4, "High", NA))
cat("Grade counts:\n"); print(table(cd$tumor_grade, useNA = "ifany"))

bt <- assay(se[unique(probe_gene$probe), tum])

# ── 3. Per-probe Spearman with grade ─────────────────────
g <- cd$grade_num
res <- do.call(rbind, lapply(rownames(bt), function(p) {
  x <- bt[p, ]; ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 30 || sd(x[ok]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], g[ok], method = "spearman", exact = FALSE))
  data.frame(probe = p, n = sum(ok), mean_beta = mean(x[ok]),
             rho_grade = unname(ct$estimate), p = ct$p.value,
             delta_high_low = mean(x[cd$grade_grp %in% "High"], na.rm = TRUE) -
               mean(x[cd$grade_grp %in% "Low"],  na.rm = TRUE))
}))
res$padj <- p.adjust(res$p, "BH")
res <- inner_join(probe_gene, res, by = "probe") %>% arrange(padj)
write.csv(res, file.path(out_dir, "probe_meth_vs_grade.csv"), row.names = FALSE)
cat("\nProbes tested:", nrow(res), "| significant (padj < 0.05):", sum(res$padj < 0.05), "\n")
cat("Significant AND |rho| >= 0.2:", sum(res$padj < 0.05 & abs(res$rho_grade) >= 0.2), "\n")
print(as.data.frame(head(res, 15)), digits = 3)
print(as.data.frame(res %>% group_by(gene) %>%
                      summarise(n_probes = n(), n_sig = sum(padj < 0.05),
                                best_rho = rho_grade[which.min(p)], best_padj = min(padj))), digits = 3)

# ── 4. Gene-level mean beta vs grade ─────────────────────
gb <- t(sapply(unique(probe_gene$gene), function(gn)
  colMeans(bt[probe_gene$probe[probe_gene$gene == gn], , drop = FALSE], na.rm = TRUE)))
gres <- do.call(rbind, lapply(rownames(gb), function(gn) {
  ok <- is.finite(gb[gn, ]) & !is.na(g)
  ct <- suppressWarnings(cor.test(gb[gn, ok], g[ok], method = "spearman", exact = FALSE))
  data.frame(gene = gn, n = sum(ok), rho_grade = unname(ct$estimate), p = ct$p.value,
             beta_G1 = mean(gb[gn, which(g == 1)], na.rm = TRUE), beta_G2 = mean(gb[gn, which(g == 2)], na.rm = TRUE),
             beta_G3 = mean(gb[gn, which(g == 3)], na.rm = TRUE), beta_G4 = mean(gb[gn, which(g == 4)], na.rm = TRUE))
}))
gres$padj <- p.adjust(gres$p, "BH")
gres <- gres %>% arrange(padj)
write.csv(gres, file.path(out_dir, "gene_meth_vs_grade.csv"), row.names = FALSE)
print(as.data.frame(gres), digits = 3)

# ── 5. Figures ───────────────────────────────────────────
res$sig <- ifelse(res$padj < 0.05, "padj < 0.05", "n.s.")
fa <- ggplot(res, aes(gene, rho_grade, color = sig)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.8) +
  scale_color_manual(values = c(`padj < 0.05` = "#D7301F", `n.s.` = "grey65"), name = NULL) +
  labs(title = "Promoter probe methylation vs grade (Spearman, G1-G4)",
       subtitle = "Each dot = one promoter probe; positive rho = more methylation in higher grade",
       x = NULL, y = "Spearman rho with grade") +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(fig_dir, "fig_probe_meth_vs_grade_rho.png"), fa, width = 10, height = 5, dpi = 300, bg = "white")

top <- head(res, 6)
d <- as.data.frame(t(bt[top$probe, , drop = FALSE])); d$grade <- cd$tumor_grade
d <- d %>% filter(grade %in% c("G1", "G2", "G3", "G4")) %>%
  pivot_longer(-grade, names_to = "probe", values_to = "beta") %>%
  left_join(top %>% dplyr::select(probe, gene), by = "probe") %>% mutate(panel = paste(gene, probe))
fb <- ggplot(d, aes(grade, beta, fill = grade)) +
  geom_boxplot(outlier.size = 0.7, show.legend = FALSE) + facet_wrap(~ panel, scales = "free_y") +
  labs(title = "Top 6 probes associated with grade", x = NULL, y = "Probe beta") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_top_probes_by_grade.png"), fb, width = 9, height = 6, dpi = 300, bg = "white")


#-----
sel <- read.csv("03_results/methylation/probes_rho_beyond_0.3.csv")
sig <- res %>% filter(padj < 0.05)
intersect(sig$probe, sel$probe)

# do the grade-associated probes also correlate with expression?
ex <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
both <- sig %>% inner_join(ex %>% select(probe, rho_expr = rho, padj_expr = padj), by = "probe") %>%
  filter(padj_expr < 0.05) %>% arrange(desc(abs(rho_expr)))
as.data.frame(both %>% select(gene, probe, rho_grade, rho_expr, padj, padj_expr, delta_high_low))