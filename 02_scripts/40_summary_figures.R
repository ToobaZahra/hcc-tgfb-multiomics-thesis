## 40: summary figures from the existing result CSVs (no raw data needed)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
m_dir <- "03_results/methylation"; d_dir <- "03_results/dge"; fig_dir <- "04_figures"

probes <- read.csv(file.path(m_dir, "probes_rho_beyond_0.3_methylation.csv"))
st     <- read.csv(file.path(m_dir, "probe_expr_cor_by_stage.csv"))
gr     <- read.csv(file.path(m_dir, "probe_expr_cor_by_grade.csv"))
gmeth  <- read.csv(file.path(m_dir, "gene_meth_by_stage_grade.csv"))
dge    <- read.csv(file.path(d_dir,  "nexus_stage_grade_dge.csv"))

col3 <- c(Hyper = "#D7301F", Hypo = "#2C7BB6", `No change` = "grey60")

# ── Fig A: probe beta, normal -> tumor ───────────────────
pa <- probes %>% mutate(label = paste(gene, probe)) %>% arrange(delta) %>%
  mutate(label = factor(label, unique(label)))
figA <- ggplot(pa, aes(y = label)) +
  geom_segment(aes(x = beta_normal, xend = beta_tumor, yend = label, color = tumor_vs_normal),
               linewidth = 0.8) +
  geom_point(aes(x = beta_normal), shape = 1, size = 2.2, color = "black") +
  geom_point(aes(x = beta_tumor, color = tumor_vs_normal), size = 2.8) +
  scale_color_manual(values = col3, name = "Tumor vs normal") +
  labs(title = "Selected probes (|rho| > 0.3): methylation in normal vs tumor",
       subtitle = "Open circle = normal mean beta; filled circle = tumor mean beta",
       x = "Mean beta", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_probe_tumor_vs_normal.png"), figA, width = 8, height = 8, dpi = 300, bg = "white")

# ── Fig B: probe rho in early vs late, low vs high ───────
b1 <- st %>% transmute(gene, probe, x = rho_early, y = rho_late, panel = "Stage: Early (x) vs Late (y)")
b2 <- gr %>% transmute(gene, probe, x = rho_gradeLow, y = rho_gradeHigh, panel = "Grade: Low (x) vs High (y)")
pb <- bind_rows(b1, b2)
figB <- ggplot(pb, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_point(size = 2.5, alpha = 0.8, color = "#2C7BB6") +
  facet_wrap(~ panel) + coord_equal() +
  labs(title = "Methylation-expression correlation of each probe inside clinical groups",
       subtitle = "Same direction in both groups; no probe differed significantly after FDR (Late/High groups are smaller, so rho is noisier)",
       x = "Spearman rho (first group)", y = "Spearman rho (second group)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_probe_rho_stage_grade.png"), figB, width = 9, height = 5, dpi = 300, bg = "white")

# ── Fig C: expression log2FC by stage and grade ──────────
pc <- dge %>% select(gene, log2FC_stage, padj_stage, log2FC_grade, padj_grade) %>%
  pivot_longer(-gene, names_to = c(".value", "comparison"), names_pattern = "(log2FC|padj)_(stage|grade)") %>%
  mutate(comparison = recode(comparison, stage = "Late vs Early stage", grade = "High vs Low grade"),
         sig = ifelse(!is.na(padj) & padj < 0.05, "padj < 0.05", "n.s."))
ord <- dge %>% arrange(log2FC_stage) %>% pull(gene)
pc$gene <- factor(pc$gene, ord)
figC <- ggplot(pc, aes(gene, log2FC, fill = sig)) +
  geom_col() + coord_flip() + facet_wrap(~ comparison) +
  scale_fill_manual(values = c(`padj < 0.05` = "#D7301F", `n.s.` = "grey70"), name = NULL) +
  labs(title = "Nexus gene expression by stage and grade (DESeq2)", x = NULL, y = "log2 fold change") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_expr_stage_grade_log2fc.png"), figC, width = 9, height = 7, dpi = 300, bg = "white")

# ── Fig D: gene-level methylation change by stage and grade ──
pd <- gmeth %>% select(gene, delta_stage, padj_stage, delta_grade, padj_grade) %>%
  pivot_longer(-gene, names_to = c(".value", "comparison"), names_pattern = "(delta|padj)_(stage|grade)") %>%
  mutate(comparison = recode(comparison, stage = "Late - Early stage", grade = "High - Low grade"),
         sig = ifelse(padj < 0.05, "padj < 0.05", "n.s."))
ordg <- gmeth %>% arrange(delta_grade) %>% pull(gene)
pd$gene <- factor(pd$gene, ordg)
figD <- ggplot(pd, aes(gene, delta, fill = sig)) +
  geom_col() + coord_flip() + facet_wrap(~ comparison) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-0.12, 0.12)) +
  scale_fill_manual(values = c(`padj < 0.05` = "#D7301F", `n.s.` = "grey70"), name = NULL) +
  labs(title = "Promoter methylation change by stage and grade (26 genes)",
       subtitle = "Dashed lines = hyper/hypo cutoff (+/-0.10); all changes are far below it",
       x = NULL, y = "Delta mean promoter beta") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_gene_meth_stage_grade_delta.png"), figD, width = 9, height = 7, dpi = 300, bg = "white")

cat("Saved 4 figures in", fig_dir, "\n")