## 37: selected probes (|rho| > 0.3) vs clinical features: early vs late stage, low vs high grade
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr); library(ggplot2); library(tidyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/methylation"; fig_dir <- "04_figures"

# ── 1. Inputs ────────────────────────────────────────────
sel <- read.csv(file.path(out_dir, "probes_rho_beyond_0.3.csv"))   # gene, probe, rho
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")

# one primary tumor per patient
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cd  <- as.data.frame(colData(se))[tum, ]
cd$patient <- substr(colnames(se)[tum], 1, 12)

# ── 2. Clinical groups ───────────────────────────────────
# Early = Stage I, II (incl. IA/IB/IIB); Late = Stage III, IV (incl. IIIA-C, IVA). Stage 0 / missing = excluded
cd$stage <- ifelse(grepl("^Stage I{1,2}[ABC]?$", cd$ajcc_pathologic_stage), "Early",
                   ifelse(grepl("^Stage (III|IV)", cd$ajcc_pathologic_stage), "Late", NA))
# Grade: G1/G2 = Low, G3/G4 = High (GX / missing excluded)
cd$grade <- ifelse(cd$tumor_grade %in% c("G1", "G2"), "Low",
                   ifelse(cd$tumor_grade %in% c("G3", "G4"), "High", NA))
cat("Tumor_grade values:\n"); print(table(cd$tumor_grade, useNA = "ifany"))
cat("\nStage groups:\n");     print(table(cd$stage, useNA = "ifany"))
cat("\nGrade groups:\n");     print(table(cd$grade, useNA = "ifany"))

# ── 3. Probe beta (tumors) ───────────────────────────────
bt <- assay(se[unique(sel$probe), tum]); colnames(bt) <- cd$patient

# ── 4. Gene expression (VST, tumors) ─────────────────────
genes <- unique(sel$gene)
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
idx <- match(sym2ens, sub("\\..*$", "", rownames(expr))); ok <- !is.na(idx)
ex <- expr[idx[ok], , drop = FALSE]; rownames(ex) <- names(sym2ens)[ok]
et <- which(substr(colnames(ex), 14, 15) == "01")
et <- et[!duplicated(substr(colnames(ex)[et], 1, 12))]
ex <- ex[, et, drop = FALSE]; colnames(ex) <- substr(colnames(ex), 1, 12)
st_ex <- setNames(cd$stage, cd$patient)[colnames(ex)]
gr_ex <- setNames(cd$grade, cd$patient)[colnames(ex)]

# ── 5. Helper: mean of two groups + Wilcoxon p ───────────
cmp <- function(v, grp, lv) {
  a <- v[which(grp == lv[1])]; b <- v[which(grp == lv[2])]
  data.frame(m1 = mean(a, na.rm = TRUE), m2 = mean(b, na.rm = TRUE),
             p = if (sum(is.finite(a)) >= 3 && sum(is.finite(b)) >= 3)
               suppressWarnings(wilcox.test(a, b)$p.value) else NA_real_)
}

# ── 6. Probe methylation vs stage and grade ──────────────
res <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  p <- sel$probe[i]
  s <- cmp(bt[p, ], cd$stage, c("Early", "Late"))
  g <- cmp(bt[p, ], cd$grade, c("Low", "High"))
  data.frame(gene = sel$gene[i], probe = p, rho = sel$rho[i],
             beta_early = s$m1, beta_late = s$m2, delta_stage = s$m2 - s$m1, p_stage = s$p,
             beta_gradeLow = g$m1, beta_gradeHigh = g$m2, delta_grade = g$m2 - g$m1, p_grade = g$p)
}))
res$padj_stage <- p.adjust(res$p_stage, "BH")
res$padj_grade <- p.adjust(res$p_grade, "BH")

# ── 7. Gene expression vs stage and grade ────────────────
resg <- do.call(rbind, lapply(intersect(genes, rownames(ex)), function(g) {
  s <- cmp(ex[g, ], st_ex, c("Early", "Late"))
  r <- cmp(ex[g, ], gr_ex, c("Low", "High"))
  data.frame(gene = g, expr_early = s$m1, expr_late = s$m2, p_expr_stage = s$p,
             expr_gradeLow = r$m1, expr_gradeHigh = r$m2, p_expr_grade = r$p)
}))
resg$padj_expr_stage <- p.adjust(resg$p_expr_stage, "BH")
resg$padj_expr_grade <- p.adjust(resg$p_expr_grade, "BH")

final <- left_join(res, resg, by = "gene") %>% arrange(rho)

# DESeq2 late vs early result (from script 09)
dge <- readRDS("03_results/dge/dge_late_vs_early.rds")
deseq <- dge %>% filter(symbol %in% genes) %>%
  group_by(symbol) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(gene = symbol, log2FC_late_vs_early = log2FoldChange, padj_deseq_stage = padj)
final <- left_join(final, deseq, by = "gene")
write.csv(final, file.path(out_dir, "probes_stage_grade.csv"), row.names = FALSE)

cat("\n== Probe methylation: early vs late stage / low vs high grade ==\n")
print(as.data.frame(final %>% dplyr::select(gene, probe, rho, delta_stage, padj_stage, delta_grade, padj_grade)))
cat("\n== Gene expression: early vs late stage / low vs high grade ==\n")
print(as.data.frame(resg %>% dplyr::select(gene, expr_early, expr_late, padj_expr_stage, expr_gradeLow, expr_gradeHigh, padj_expr_grade)))

# ── 8. Boxplots for probes with nominal p < 0.05 (stage) ─
hit <- final %>% filter(p_stage < 0.05) %>% pull(probe)
if (length(hit) > 0) {
  lab <- final %>% dplyr::select(probe, gene) %>% distinct()
  d <- as.data.frame(t(bt[hit, , drop = FALSE])); d$patient <- rownames(d)
  d$stage <- cd$stage[match(d$patient, cd$patient)]
  d <- d %>% filter(!is.na(stage)) %>% pivot_longer(all_of(hit), names_to = "probe", values_to = "beta") %>%
    left_join(lab, by = "probe") %>% mutate(panel = paste(gene, probe))
  p <- ggplot(d, aes(factor(stage, c("Early", "Late")), beta, fill = stage)) +
    geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
    facet_wrap(~ panel, scales = "free_y") +
    labs(x = NULL, y = "Probe beta", title = "Probes with stage difference (p < 0.05)") +
    theme_minimal(base_size = 11)
  ggsave(file.path(fig_dir, "fig_probe_stage_boxplot.png"), p, width = 10, height = 7, dpi = 300, bg = "white")
}

d <- as.data.frame(t(ex)); d$patient <- rownames(d)
d$Stage <- st_ex; d$Grade <- gr_ex
long <- d %>% pivot_longer(-c(patient, Stage, Grade), names_to = "gene", values_to = "expr")

# panel labels with padj from resg
lab_g <- resg %>% mutate(label = paste0(gene, "\npadj = ", signif(padj_expr_grade, 2))) %>%
  dplyr::select(gene, label)
lab_s <- resg %>% mutate(label = paste0(gene, "\npadj = ", signif(padj_expr_stage, 2))) %>%
  dplyr::select(gene, label)

long_g <- long %>% filter(!is.na(Grade)) %>% left_join(lab_g, by = "gene")
long_s <- long %>% filter(!is.na(Stage)) %>% left_join(lab_s, by = "gene")

p_grade <- ggplot(long_g, aes(factor(Grade, c("Low", "High")), expr, fill = Grade)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = "Grade", y = "VST expression", title = "Expression by grade (Low = G1/G2, High = G3/G4)") +
  theme_minimal(base_size = 11)

p_stage <- ggplot(long_s, aes(factor(Stage, c("Early", "Late")), expr, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  facet_wrap(~ label, scales = "free_y") +
  labs(x = "Stage", y = "VST expression", title = "Expression by stage (Early = I/II, Late = III/IV)") +
  theme_minimal(base_size = 11)

ggsave("04_figures/fig_expr_by_grade.png", p_grade, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("04_figures/fig_expr_by_stage.png", p_stage, width = 10, height = 8, dpi = 300, bg = "white")



#-----------------------------------------------------
pts <- intersect(colnames(bt), colnames(ex))
grp_of <- function(col) setNames(cd[[col]], cd$patient)[pts]
rho_in <- function(p, g, k) {
  x <- bt[p, k]; y <- ex[g, k]; ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(NA, NA, sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(unname(ct$estimate), ct$p.value, sum(ok))
}
by_group <- function(col, lv, tag) {
  gp <- grp_of(col)
  do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
    if (!sel$gene[i] %in% rownames(ex)) return(NULL)
    a <- rho_in(sel$probe[i], sel$gene[i], pts[which(gp == lv[1])])
    b <- rho_in(sel$probe[i], sel$gene[i], pts[which(gp == lv[2])])
    z <- (atanh(a[1]) - atanh(b[1])) / sqrt(1/(a[3]-3) + 1/(b[3]-3))
    out <- data.frame(gene = sel$gene[i], probe = sel$probe[i], rho_all = sel$rho[i],
                      r1 = a[1], p1 = a[2], n1 = a[3], r2 = b[1], p2 = b[2], n2 = b[3],
                      p_diff = 2 * pnorm(-abs(z)))
    names(out)[4:9] <- paste0(c("rho_", "p_", "n_"), rep(tag, each = 3))
    out
  }))
}
st_cor <- by_group("stage", c("Early", "Late"), c("early", "late"))
gr_cor <- by_group("grade", c("Low", "High"), c("gradeLow", "gradeHigh"))
st_cor$padj_diff <- p.adjust(st_cor$p_diff, "BH")
gr_cor$padj_diff <- p.adjust(gr_cor$p_diff, "BH")
write.csv(st_cor, file.path(out_dir, "probe_expr_cor_by_stage.csv"), row.names = FALSE)
write.csv(gr_cor, file.path(out_dir, "probe_expr_cor_by_grade.csv"), row.names = FALSE)
print(as.data.frame(st_cor))