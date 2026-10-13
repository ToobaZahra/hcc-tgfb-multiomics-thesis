## 47: final integrated table -- one row per gene combining expression, methylation, grade, stage, survival
## (reads existing result files only; run scripts 43, 44, 45, 46 first)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; m <- "03_results/methylation"; fig_dir <- "04_figures"
sig <- function(p) !is.na(p) & p < 0.05

# ── 1. Base: summary from script 44 (grade, stage, DESeq2, tumor-normal, link) ──
S <- read.csv(file.path(a, "7_summary_best_findings.csv"))

# ── 2. Expression vs survival (script 45) ────────────────
es <- read.csv(file.path(a, "8_expr_vs_survival_gene.csv")) %>%
  transmute(gene, expr_HR = HR, expr_surv_padj = padj_cox,
            expr_HR_adj = HR_adj, expr_surv_p_adj_model = p_adj_model, expr_logrank_padj = padj_logrank)

# ── 3. Methylation vs survival, gene level (script 46) ───
ms <- read.csv(file.path(a, "9_meth_vs_survival_gene.csv")) %>%
  transmute(gene, meth_HR = HR, meth_surv_padj = padj_cox,
            meth_HR_adj = HR_adj, meth_surv_p_adj_model = p_adj_model)

# ── 4. Methylation vs survival, probe level ──────────────
tp <- read.csv(file.path(a, "10_meth_vs_survival_probe.csv"))
pg <- read.csv(file.path(a, "5_meth_vs_grade_probe.csv")) %>% filter(significant) %>% dplyr::select(gene, probe)
pe <- read.csv(file.path(m, "meth_expr_corr_probe_level.csv")) %>% filter(sig(padj)) %>% dplyr::select(gene, probe)
ts <- tp %>% filter(sig(padj_cox))
ps <- tp %>% group_by(gene) %>%
  summarise(n_probes_surv_sig = sum(sig(padj_cox)), best_probe_surv = probe[which.min(p_cox)],
            best_probe_HR = HR[which.min(p_cox)], .groups = "drop")
b1 <- ts %>% inner_join(pe, by = c("gene", "probe")) %>% count(gene, name = "n_probes_surv_and_expr")
b2 <- ts %>% inner_join(pg, by = c("gene", "probe")) %>% count(gene, name = "n_probes_surv_and_grade")

F <- S %>% left_join(es, by = "gene") %>% left_join(ms, by = "gene") %>%
  left_join(ps, by = "gene") %>% left_join(b1, by = "gene") %>% left_join(b2, by = "gene")
for (v in c("n_probes_surv_sig", "n_probes_surv_and_expr", "n_probes_surv_and_grade"))
  F[[v]][is.na(F[[v]])] <- 0

# ── 5. Evidence flags (padj < 0.05) and tiers ────────────
flags <- data.frame(
  `Expr~grade`      = sig(F$expr_padj_grade),
  `Expr~stage`      = sig(F$expr_padj_stage),
  `Meth~grade`      = sig(F$meth_padj_grade),
  `Meth~stage`      = sig(F$meth_padj_stage),
  `DESeq2 grade`    = sig(F$deseq_padj_grade),
  `DESeq2 stage`    = sig(F$deseq_padj_stage),
  `Expr~survival`   = sig(F$expr_surv_padj),
  `Meth~survival`   = sig(F$meth_surv_padj),
  `Probe~survival`  = F$n_probes_surv_sig > 0,
  check.names = FALSE)
F$n_evidence_total <- rowSums(flags)
F$significant_in <- apply(flags, 1, function(r) paste(names(flags)[r], collapse = "; "))
F$tier <- case_when(F$n_evidence_total >= 4 ~ "Tier 1 (4+ analyses)",
                    F$n_evidence_total >= 2 ~ "Tier 2 (2-3 analyses)",
                    TRUE ~ "Tier 3 (0-1 analyses)")
F <- F %>% arrange(desc(n_evidence_total), expr_surv_padj)
write.csv(F, file.path(a, "11_final_summary.csv"), row.names = FALSE)

cat("\n== Final summary (max 9 analyses) ==\n")
print(as.data.frame(F %>% dplyr::select(gene, n_evidence_total, tier, significant_in)), right = FALSE)
cat("\nKey numbers per gene:\n")
print(as.data.frame(F %>% dplyr::select(gene, expr_rho_grade, expr_rho_stage, meth_rho_grade, expr_HR, meth_HR,
                                        n_probes_surv_sig, n_probes_surv_and_expr, n_probes_surv_and_grade) %>% head(12)), digits = 2)

# ── 6. Figure: which analyses are significant for which gene (colour = direction) ──
dirs <- data.frame(
  gene = F$gene,
  `Expr~grade` = sign(F$expr_rho_grade), `Expr~stage` = sign(F$expr_rho_stage),
  `Meth~grade` = sign(F$meth_rho_grade), `Meth~stage` = sign(F$meth_rho_stage),
  `DESeq2 grade` = sign(F$deseq_log2FC_grade), `DESeq2 stage` = sign(F$deseq_log2FC_stage),
  `Expr~survival` = sign(log(F$expr_HR)), `Meth~survival` = sign(log(F$meth_HR)),
  `Probe~survival` = sign(log(F$best_probe_HR)), check.names = FALSE)
long <- dirs %>% pivot_longer(-gene, names_to = "analysis", values_to = "dir") %>%
  left_join(flags %>% mutate(gene = F$gene) %>% pivot_longer(-gene, names_to = "analysis", values_to = "is_sig"),
            by = c("gene", "analysis")) %>%
  mutate(cell = case_when(!is_sig ~ "n.s.", dir > 0 ~ "Positive / higher", TRUE ~ "Negative / lower"))
long$gene <- factor(long$gene, rev(F$gene)); long$analysis <- factor(long$analysis, names(flags))
f <- ggplot(long, aes(analysis, gene, fill = cell)) + geom_tile(color = "white") +
  scale_fill_manual(values = c(`n.s.` = "grey92", `Positive / higher` = "#D7301F", `Negative / lower` = "#2C7BB6"), name = NULL) +
  labs(title = "Where each nexus gene shows a significant association (padj < 0.05)",
       subtitle = "Red = higher expression/methylation/hazard with grade, stage or risk; blue = lower",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "fig_final_summary.png"), f, width = 8, height = 8, dpi = 300, bg = "white")