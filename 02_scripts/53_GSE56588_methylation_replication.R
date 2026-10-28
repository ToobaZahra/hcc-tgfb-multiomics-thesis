## 53: DUSP10 cross-cohort summary -- TCGA-LIHC + GSE14520 + GSE76427 + GSE56588
## Combines outputs of script 51 (files 20, 21, 22) into one long and one wide table.
## Expression: TCGA, GSE14520, GSE76427.  Methylation: TCGA, GSE56588.
## GSE56588 rows are added from the GSE56588 script outputs (files 17, 17b) if they exist.
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; fig_dir <- "04_figures"
t20 <- read.csv(file.path(a, "20_DUSP10_survival_models.csv"))
t21 <- read.csv(file.path(a, "21_DUSP10_TCGA_tests.csv"))
t22 <- read.csv(file.path(a, "22_DUSP10_GEO_tests.csv"))

# one row of the combined table
mk <- function(layer, analysis, dataset, estimate, p, lo = NA, hi = NA, stat, n = NA) {
  data.frame(layer = layer, analysis = analysis, dataset = dataset, statistic = stat,
             estimate = estimate, ci_low = lo, ci_high = hi, p = p, n = n,
             stringsAsFactors = FALSE)
}

# ---- TCGA-LIHC (from files 20 and 21) ----
g21 <- function(f, t) t21[t21$feature == f & t21$test == t, ]
g20 <- function(m, term) t20[t20$model == m & t20$term == term, ]
tc <- "TCGA-LIHC"

x <- g21("Expression", "Tumor minus normal (mean VST)")
r_tn <- mk("Expression", "Tumor vs normal", tc, x$estimate, x$p, stat = "mean difference", n = x$n)
x <- g21("Expression", "Late minus Early stage (median VST)")
r_le <- mk("Expression", "Late vs early stage", tc, x$estimate, x$p, stat = "median difference", n = x$n)
x <- g21("Expression", "Stage I-IV ordered (Spearman rho)")
r_st <- mk("Expression", "Stage I-IV trend", tc, x$estimate, x$p, stat = "Spearman rho", n = x$n)
x <- g21("Expression", "High minus Low grade (median VST)")
r_hl <- mk("Expression", "High vs low grade", tc, x$estimate, x$p, stat = "median difference", n = x$n)
x <- g21("Expression", "Grade G1-G4 ordered (Spearman rho)")
r_gr <- mk("Expression", "Grade G1-G4 trend", tc, x$estimate, x$p, stat = "Spearman rho", n = x$n)
x <- g20("Expression alone", "expr_z")
r_hr <- mk("Expression", "Survival, unadjusted", tc, x$HR, x$p, x$lo, x$hi, stat = "HR per 1 SD", n = x$n)
x <- g20("Expression + age + stage", "expr_z")
r_ha <- mk("Expression", "Survival, adjusted (age + stage)", tc, x$HR, x$p, x$lo, x$hi, stat = "HR per 1 SD", n = x$n)

x <- g21("Gene methylation", "Tumor minus normal (delta beta)")
m_tn <- mk("Methylation", "Tumor vs normal", tc, x$estimate, x$p, stat = "delta beta", n = x$n)
x <- g21("Gene methylation", "Late minus Early stage (median beta)")
m_le <- mk("Methylation", "Late vs early stage", tc, x$estimate, x$p, stat = "median difference", n = x$n)
x <- g21("Gene methylation", "Grade G1-G4 ordered (Spearman rho)")
m_gr <- mk("Methylation", "Grade G1-G4 trend", tc, x$estimate, x$p, stat = "Spearman rho", n = x$n)
x <- g21("Methylation vs expression", "All tumors (Spearman rho)")
m_ex <- mk("Methylation", "Methylation vs expression", tc, x$estimate, x$p, stat = "Spearman rho", n = x$n)
x <- g20("Methylation alone", "meth_z")
m_hr <- mk("Methylation", "Survival, unadjusted", tc, x$HR, x$p, x$lo, x$hi, stat = "HR per 1 SD", n = x$n)

tcga <- bind_rows(r_tn, r_le, r_st, r_hl, r_gr, r_hr, r_ha, m_tn, m_le, m_gr, m_ex, m_hr)

# ---- GEO cohorts (from file 22, expression only) ----
geo_rows <- function(r) bind_rows(
  mk("Expression", "Tumor vs normal", r$dataset, r$tn_diff, r$tn_p, stat = "mean difference"),
  mk("Expression", "Late vs early stage", r$dataset, r$late_early_diff, r$late_early_p, stat = "median difference"),
  mk("Expression", "Stage I-IV trend", r$dataset, r$stage_rho, r$stage_p, stat = "Spearman rho"),
  mk("Expression", "Survival, unadjusted", r$dataset, r$HR, r$cox_p, r$HR_low, r$HR_high,
     stat = "HR per 1 SD", n = r$n_surv),
  mk("Expression", "Survival, adjusted (age + stage)", r$dataset, r$HR_adj, r$cox_p_adj,
     stat = "HR per 1 SD", n = r$n_surv))
geo <- bind_rows(lapply(split(t22, seq_len(nrow(t22))), geo_rows))

# ---- GSE56588 (methylation only; from the GSE56588 script outputs) ----
f17 <- file.path(a, "17_GSE56588_meth_7genes.csv")
f17b <- file.path(a, "17b_GSE56588_probe_delta.csv")
gse56588 <- NULL
if (file.exists(f17)) {
  t17 <- read.csv(f17) %>% filter(gene == "DUSP10")
  gse56588 <- mk("Methylation", "Tumor vs normal", "GSE56588", t17$delta_GSE56588, t17$p_wilcox,
                 stat = "delta beta")
  cat("GSE56588 DUSP10: delta =", round(t17$delta_GSE56588, 4), "| Holm p over 7 genes =",
      signif(t17$p_holm, 3), "\n")
} else cat("File 17 not found: run the GSE56588 script first. GSE56588 column will be empty.\n")

# ---- combine ----
long <- bind_rows(tcga, geo, gse56588) %>%
  mutate(dataset = factor(dataset, c("TCGA-LIHC", "GSE14520", "GSE76427", "GSE56588")),
         direction = ifelse(grepl("^HR", statistic), ifelse(estimate < 1, "lower risk", "higher risk"),
                            ifelse(estimate < 0, "down / negative", "up / positive")),
         sig_raw = !is.na(p) & p < 0.05) %>%
  arrange(layer, analysis, dataset)
write.csv(long, file.path(a, "23_DUSP10_cross_cohort_long.csv"), row.names = FALSE)

wide <- long %>%
  mutate(txt = ifelse(is.na(ci_low),
                      sprintf("%.3g (p=%.2g)", estimate, p),
                      sprintf("%.2f [%.2f-%.2f] (p=%.2g)", estimate, ci_low, ci_high, p))) %>%
  select(layer, analysis, statistic, dataset, txt) %>%
  pivot_wider(names_from = dataset, values_from = txt)
write.csv(wide, file.path(a, "23_DUSP10_cross_cohort_summary.csv"), row.names = FALSE)

print(as.data.frame(wide), right = FALSE)
cat("\nSaved 23_DUSP10_cross_cohort_summary.csv (readable) and 23_DUSP10_cross_cohort_long.csv (for plotting).\n")

# ---- probe-level methylation: TCGA vs GSE56588 ----
f19b <- file.path(a, "19b_DUSP10_probe_tumor_vs_normal.csv")
if (file.exists(f17b) && file.exists(f19b)) {
  p19b <- read.csv(f19b) %>% select(probe, delta_TCGA = delta)
  p17b <- read.csv(f17b) %>% filter(gene == "DUSP10") %>% select(probe, delta_GSE56588 = delta)
  p19 <- read.csv(file.path(a, "19_DUSP10_probe_table.csv")) %>% select(probe, dist_TSS, rho_expr, padj_expr)
  pl <- p19b %>% left_join(p17b, by = "probe") %>% left_join(p19, by = "probe") %>%
    mutate(same_direction = sign(delta_TCGA) == sign(delta_GSE56588)) %>% arrange(dist_TSS)
  write.csv(pl, file.path(a, "23b_DUSP10_methylation_probe_level.csv"), row.names = FALSE)
  cat("\n== DUSP10 probe-level tumor minus normal (TCGA vs GSE56588) ==\n"); print(pl, digits = 3)
  cat("Probes missing in GSE56588:", sum(is.na(pl$delta_GSE56588)), "of", nrow(pl),
      "| same direction:", sum(pl$same_direction, na.rm = TRUE), "of", sum(!is.na(pl$delta_GSE56588)), "\n")
  fp <- ggplot(pl[!is.na(pl$delta_GSE56588), ], aes(delta_TCGA, delta_GSE56588, label = probe)) +
    geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = dist_TSS), size = 3) + geom_text(size = 2.8, vjust = -0.9, check_overlap = TRUE) +
    scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7301F", midpoint = 0, name = "Distance\nto TSS (bp)") +
    labs(title = "DUSP10 promoter probes: tumor minus normal methylation",
         subtitle = "Each dot = one CpG probe; dashed line = identical delta in both cohorts",
         x = "Delta beta, TCGA-LIHC", y = "Delta beta, GSE56588") + theme_minimal(base_size = 11)
  ggsave(file.path(fig_dir, "fig_DUSP10_meth_TCGA_vs_GSE56588.png"), fp, width = 6.5, height = 5.5, dpi = 300, bg = "white")
  cat("Saved 23b_DUSP10_methylation_probe_level.csv and fig_DUSP10_meth_TCGA_vs_GSE56588.png\n")
} else cat("\nProbe-level section skipped: file 17b or 19b not found.\n")

cat("\nNotes: p-values are raw. TCGA expression is VST, GEO is log2, so compare direction, not magnitude.\n")
cat("GSE56588 has ~10 normals, so rely on delta and direction more than p. GSE56588 holds methylation only;\n")
cat("methylation vs expression, grade, stage and survival remain TCGA-only. Expression rows are NA for GSE56588.\n")