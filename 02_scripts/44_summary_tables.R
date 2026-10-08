## 44: one summary table of the best findings per gene (reads existing result files only)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; m <- "03_results/methylation"; d <- "03_results/dge"; fig_dir <- "04_figures"
sig <- function(p) !is.na(p) & p < 0.05

# ── 1. Gene-level association tables (script 43) ─────────
t1 <- read.csv(file.path(a, "1_expr_vs_grade_gene.csv"))
t2 <- read.csv(file.path(a, "2_expr_vs_stage_gene.csv"))
t3 <- read.csv(file.path(a, "3_meth_vs_grade_gene.csv"))
t4 <- read.csv(file.path(a, "4_meth_vs_stage_gene.csv"))
eg <- t1 %>% transmute(gene, expr_rho_grade = rho, expr_padj_grade = padj, expr_G1_to_G4 = expr_G4 - expr_G1)
es <- t2 %>% transmute(gene, expr_rho_stage = rho, expr_padj_stage = padj)
mg <- t3 %>% transmute(gene, meth_rho_grade = rho, meth_padj_grade = padj, beta_G1_to_G4 = beta_G4 - beta_G1)
ms <- t4 %>% transmute(gene, meth_rho_stage = rho, meth_padj_stage = padj)

# ── 2. Probe-level summary per gene ──────────────────────
pg <- read.csv(file.path(a, "5_meth_vs_grade_probe.csv"))
ps <- read.csv(file.path(a, "6_meth_vs_stage_probe.csv"))
psum <- function(x, tag) {
  o <- x %>% group_by(gene) %>%
    summarise(n_probes = n(), n_sig = sum(significant),
              best_probe = probe[which.min(p)], best_rho = rho[which.min(p)], .groups = "drop")
  names(o)[-1] <- paste0(names(o)[-1], "_", tag); o
}
pe <- read.csv(file.path(m, "meth_expr_corr_probe_level.csv")) %>% filter(sig(padj)) %>% dplyr::select(gene, probe)
bg <- pg %>% filter(significant) %>% inner_join(pe, by = c("gene", "probe")) %>% count(gene, name = "n_probes_grade_and_expr")
bs <- ps %>% filter(significant) %>% inner_join(pe, by = c("gene", "probe")) %>% count(gene, name = "n_probes_stage_and_expr")

# ── 3. DESeq2 and tumor-vs-normal methylation ────────────
dg <- read.csv(file.path(d, "nexus_stage_grade_dge.csv")) %>%
  transmute(gene, deseq_log2FC_tumor_normal = log2FC_tumor_vs_normal,
            deseq_log2FC_stage = log2FC_stage, deseq_padj_stage = padj_stage,
            deseq_log2FC_grade = log2FC_grade, deseq_padj_grade = padj_grade)
tn_file <- file.path(m, "meth_promoter_gene_summary.csv")
tn <- if (file.exists(tn_file)) read.csv(tn_file) %>%
  transmute(gene, meth_delta_tumor_normal = delta, pct_hyper, pct_hypo) else NULL
lk_file <- file.path(m, "meth_expr_corr_promoter_tumor.csv")
lk <- if (file.exists(lk_file)) { x <- read.csv(lk_file)
if (all(c("gene", "rho", "padj") %in% names(x))) transmute(x, gene, link_rho = rho, link_padj = padj) else NULL } else NULL

# ── 4. Join everything ───────────────────────────────────
S <- Reduce(function(x, y) left_join(x, y, by = "gene"),
            Filter(Negate(is.null), list(eg, es, mg, ms, psum(pg, "grade"), psum(ps, "stage"),
                                         bg, bs, dg, tn, lk)))
S[is.na(S$n_probes_grade_and_expr), "n_probes_grade_and_expr"] <- 0
S[is.na(S$n_probes_stage_and_expr), "n_probes_stage_and_expr"] <- 0

# ── 5. Patterns and evidence score ───────────────────────
pattern <- function(mp, ep) case_when(sig(mp) & sig(ep) ~ "Methylation + expression",
                                      sig(ep) ~ "Expression only",
                                      sig(mp) ~ "Methylation only", TRUE ~ "None")
S$grade_pattern <- pattern(S$meth_padj_grade, S$expr_padj_grade)
S$stage_pattern <- pattern(S$meth_padj_stage, S$expr_padj_stage)
S$n_evidence <- sig(S$expr_padj_grade) + sig(S$expr_padj_stage) + sig(S$meth_padj_grade) +
  sig(S$meth_padj_stage) + sig(S$deseq_padj_grade) + sig(S$deseq_padj_stage)
if ("link_rho" %in% names(S)) {
  # expected expression direction if methylation acts through the gene-level methylation-expression link
  S$grade_expected_vs_observed <- ifelse(sig(S$meth_padj_grade) & sig(S$expr_padj_grade),
                                         ifelse(sign(S$meth_rho_grade) * sign(S$link_rho) == sign(S$expr_rho_grade), "concordant", "discordant"), NA)
}
S <- S %>% arrange(desc(n_evidence), meth_padj_grade)
write.csv(S, file.path(a, "7_summary_best_findings.csv"), row.names = FALSE)

cat("\n== Summary: main columns ==\n")
print(as.data.frame(S %>% dplyr::select(gene, n_evidence, grade_pattern, stage_pattern,
                                        expr_rho_grade, meth_rho_grade, expr_rho_stage, meth_rho_stage,
                                        n_sig_grade, n_sig_stage, n_probes_grade_and_expr, best_probe_grade)), digits = 2)
cat("\nGenes with methylation + expression change with grade:",
    paste(S$gene[S$grade_pattern == "Methylation + expression"], collapse = ", "), "\n")
cat("Genes with methylation + expression change with stage:",
    paste(S$gene[S$stage_pattern == "Methylation + expression"], collapse = ", "), "\n")

# ── 6. Overview figure: rho of the four comparisons ──────
pl <- S %>% dplyr::select(gene, n_evidence,
                          `Expression ~ grade` = expr_rho_grade, `Expression ~ stage` = expr_rho_stage,
                          `Methylation ~ grade` = meth_rho_grade, `Methylation ~ stage` = meth_rho_stage) %>%
  pivot_longer(-c(gene, n_evidence), names_to = "test", values_to = "rho")
pp <- S %>% dplyr::select(gene, `Expression ~ grade` = expr_padj_grade, `Expression ~ stage` = expr_padj_stage,
                          `Methylation ~ grade` = meth_padj_grade, `Methylation ~ stage` = meth_padj_stage) %>%
  pivot_longer(-gene, names_to = "test", values_to = "padj")
pl <- left_join(pl, pp, by = c("gene", "test")) %>% mutate(star = ifelse(sig(padj), "*", ""))
pl$gene <- factor(pl$gene, rev(S$gene))
f <- ggplot(pl, aes(test, gene, fill = rho)) + geom_tile(color = "white") +
  geom_text(aes(label = star), size = 6) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7301F", midpoint = 0, name = "Spearman rho") +
  labs(title = "Association of expression and methylation with grade and stage",
       subtitle = "* = padj < 0.05; genes ordered by number of significant results", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "fig_association_summary.png"), f, width = 7, height = 8, dpi = 300, bg = "white")