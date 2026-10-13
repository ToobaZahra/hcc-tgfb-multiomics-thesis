## 48: SMAD2 deep dive -- methylation, expression, grade, stage and survival in one place
## Part 1 uses existing result CSVs only; Part 2 uses the raw TCGA-LIHC objects for SMAD2.
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(survival)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; m <- "03_results/methylation"; fig_dir <- "04_figures"
G <- "SMAD2"
sig <- function(p) !is.na(p) & p < 0.05

# ═════════ PART 1: promoter probe profile (from existing CSVs) ═════════
pe <- read.csv(file.path(m, "meth_expr_corr_probe_level.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, dist_TSS, mean_beta, rho_expr = rho, padj_expr = padj)
pg <- read.csv(file.path(a, "5_meth_vs_grade_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, rho_grade = rho, padj_grade = padj)
ps <- read.csv(file.path(a, "6_meth_vs_stage_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, rho_stage = rho, padj_stage = padj)
pv <- read.csv(file.path(a, "10_meth_vs_survival_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, HR_surv = HR, padj_surv = padj_cox)
P <- Reduce(function(x, y) left_join(x, y, by = "probe"), list(pe, pg, ps, pv)) %>% arrange(dist_TSS)
write.csv(P, file.path(a, "12_SMAD2_probe_table.csv"), row.names = FALSE)
cat("== SMAD2 promoter probes ==\n"); print(as.data.frame(P), digits = 2)

pl <- P %>% transmute(probe, dist_TSS,
                      `rho with expression` = rho_expr, `rho with grade` = rho_grade,
                      `rho with stage` = rho_stage, `log2 HR (survival)` = log2(HR_surv)) %>%
  pivot_longer(-c(probe, dist_TSS), names_to = "metric", values_to = "value")
ps_long <- P %>% transmute(probe, `rho with expression` = padj_expr, `rho with grade` = padj_grade,
                           `rho with stage` = padj_stage, `log2 HR (survival)` = padj_surv) %>%
  pivot_longer(-probe, names_to = "metric", values_to = "padj")
pl <- left_join(pl, ps_long, by = c("probe", "metric")) %>% mutate(sig = ifelse(sig(padj), "padj < 0.05", "n.s."))
fA <- ggplot(pl, aes(dist_TSS, value, color = sig)) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) + facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(`padj < 0.05` = "#D7301F", `n.s.` = "grey65"), name = NULL) +
  labs(title = "SMAD2 promoter probes: expression, grade, stage and survival",
       subtitle = "Each dot = one CpG probe; dashed line = TSS", x = "Distance to TSS (bp)", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_SMAD2_probe_profile.png"), fA, width = 7, height = 9, dpi = 300, bg = "white")

# ═════════ PART 2: patient-level analysis (raw data) ═════════
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cdm <- as.data.frame(colData(se))[tum, ]; cdm$patient <- substr(colnames(se)[tum], 1, 12)

probes <- P$probe
bt <- assay(se[probes, tum]); colnames(bt) <- cdm$patient
meth <- colMeans(bt, na.rm = TRUE)                                   # gene-level promoter beta

vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens  <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row  <- match(ens, sub("\\..*$", "", rownames(expr)))
et   <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ex   <- expr[row, et]; names(ex) <- substr(colnames(expr)[et], 1, 12)

age_col <- intersect(c("age_at_index", "age_at_diagnosis"), colnames(cdm))[1]
age <- as.numeric(cdm[[age_col]]); if (median(age, na.rm = TRUE) > 200) age <- age / 365.25
dead <- tolower(cdm$vital_status) == "dead"
time <- as.numeric(ifelse(dead, cdm$days_to_death, cdm$days_to_last_follow_up)) / 30.44
st <- cdm$ajcc_pathologic_stage
d <- data.frame(patient = cdm$patient, meth = meth[cdm$patient], expr = ex[cdm$patient],
                grade = factor(cdm$tumor_grade, c("G1", "G2", "G3", "G4")), time = time, event = as.integer(dead), age = age,
                stage = factor(ifelse(grepl("^Stage I{1,2}[ABC]?$", st), "Early",
                                      ifelse(grepl("^Stage (III|IV)", st), "Late", NA)), c("Early", "Late")))
d$grade_num <- as.integer(d$grade)
cat("\nPatients with methylation + expression:", sum(is.finite(d$meth) & is.finite(d$expr)), "\n")

# 2a. methylation vs expression: overall, adjusted for grade, within grade groups
ok <- is.finite(d$meth) & is.finite(d$expr)
cat("\n== SMAD2 methylation vs expression ==\n")
ct <- cor.test(d$meth[ok], d$expr[ok], method = "spearman", exact = FALSE)
cat(sprintf("All tumors: rho = %.3f, p = %.2g (n = %d)\n", ct$estimate, ct$p.value, sum(ok)))
k <- ok & !is.na(d$grade_num)
rm_ <- resid(lm(meth ~ grade_num, data = d[k, ])); re_ <- resid(lm(expr ~ grade_num, data = d[k, ]))
ct2 <- cor.test(rm_, re_, method = "spearman", exact = FALSE)
cat(sprintf("Adjusted for grade (residuals): rho = %.3f, p = %.2g (n = %d)\n", ct2$estimate, ct2$p.value, sum(k)))
for (gr in c("G1", "G2", "G3", "G4")) {
  kk <- ok & d$grade %in% gr
  if (sum(kk) >= 10) { c3 <- cor.test(d$meth[kk], d$expr[kk], method = "spearman", exact = FALSE)
  cat(sprintf("  %s: rho = %.3f, p = %.2g (n = %d)\n", gr, c3$estimate, c3$p.value, sum(kk))) }
}

# 2b. survival: expression, methylation, both, both + age + stage
d$expr_z <- as.numeric(scale(d$expr)); d$meth_z <- as.numeric(scale(d$meth))
ds <- d[is.finite(d$time) & d$time > 0 & is.finite(d$expr_z) & is.finite(d$meth_z), ]
cox_out <- function(f, label) {
  fit <- coxph(f, data = ds); s <- summary(fit)
  data.frame(model = label, term = rownames(s$coefficients), HR = s$conf.int[, 1], lo = s$conf.int[, 3],
             hi = s$conf.int[, 4], p = s$coefficients[, 5], n = fit$n, events = fit$nevent, row.names = NULL)
}
models <- rbind(
  cox_out(Surv(time, event) ~ expr_z, "Expression alone"),
  cox_out(Surv(time, event) ~ meth_z, "Methylation alone"),
  cox_out(Surv(time, event) ~ expr_z + meth_z, "Expression + methylation"),
  cox_out(Surv(time, event) ~ expr_z + meth_z + age + stage, "Expression + methylation + age + stage"))
write.csv(models, file.path(a, "13_SMAD2_survival_models.csv"), row.names = FALSE)
cat("\n== SMAD2 survival models (HR per 1 SD) ==\n"); print(models, digits = 3)

# 2c. figures
long <- d %>% filter(!is.na(grade)) %>% dplyr::select(grade, `Promoter methylation (beta)` = meth, `Expression (VST)` = expr) %>%
  pivot_longer(-grade, names_to = "measure", values_to = "value")
fB <- ggplot(long, aes(grade, value, fill = grade)) + geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  facet_wrap(~ measure, scales = "free_y") + labs(title = "SMAD2 by tumor grade", x = NULL, y = NULL) + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_SMAD2_by_grade.png"), fB, width = 8, height = 4, dpi = 300, bg = "white")

fC <- ggplot(d[ok, ], aes(meth, expr, color = grade)) + geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE) +
  labs(title = "SMAD2: promoter methylation vs expression", x = "Mean promoter beta", y = "VST expression") + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_SMAD2_meth_vs_expr.png"), fC, width = 6, height = 5, dpi = 300, bg = "white")

km <- function(x, nm) {
  dd <- ds; dd$grp <- factor(ifelse(x > median(x), "High", "Low"), c("Low", "High"))
  lr <- survdiff(Surv(time, event) ~ grp, data = dd)
  plot(survfit(Surv(time, event) ~ grp, data = dd), col = c("#2C7BB6", "#D7301F"), lwd = 2, xlab = "Months", ylab = "Overall survival",
       main = sprintf("%s\nlog-rank p = %.3g", nm, 1 - pchisq(lr$chisq, 1)), cex.main = 0.9)
  legend("bottomleft", c("Low", "High"), col = c("#2C7BB6", "#D7301F"), lwd = 2, bty = "n", cex = 0.8)
}
png(file.path(fig_dir, "fig_SMAD2_km.png"), width = 2000, height = 1000, res = 250)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); km(ds$expr, "SMAD2 expression"); km(ds$meth, "SMAD2 promoter methylation"); dev.off()