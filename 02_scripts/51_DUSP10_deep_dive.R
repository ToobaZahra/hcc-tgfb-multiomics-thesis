## 51: DUSP10 deep dive -- expression, methylation, stage, grade, survival (TCGA) + GEO replication
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(GEOquery)
library(survival)
library(dplyr); library(tidyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; m <- "03_results/methylation"; d_dir <- "03_results/dge"
fig_dir <- "04_figures"; g <- "data/raw/geo/"
G <- "DUSP10"
sig <- function(p) !is.na(p) & p < 0.05

# ═════════ PART 1: promoter probe profile (existing CSVs) ═════════
pe <- read.csv(file.path(m, "meth_expr_corr_probe_level.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, dist_TSS, mean_beta, rho_expr = rho, padj_expr = padj)
pg <- read.csv(file.path(a, "5_meth_vs_grade_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, rho_grade = rho, padj_grade = padj)
ps <- read.csv(file.path(a, "6_meth_vs_stage_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, rho_stage = rho, padj_stage = padj)
pv <- read.csv(file.path(a, "10_meth_vs_survival_probe.csv")) %>% filter(gene == G) %>%
  dplyr::select(probe, HR_surv = HR, padj_surv = padj_cox)
P <- Reduce(function(x, y) left_join(x, y, by = "probe"), list(pe, pg, ps, pv)) %>% arrange(dist_TSS)
write.csv(P, file.path(a, "19_DUSP10_probe_table.csv"), row.names = FALSE)
cat("== DUSP10 promoter probes ==\n"); print(as.data.frame(P), digits = 2)

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
  labs(title = "DUSP10 promoter probes: expression, grade, stage and survival",
       subtitle = "Each dot = one CpG probe; dashed line = TSS", x = "Distance to TSS (bp)", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_DUSP10_probe_profile.png"), fA, width = 7, height = 9, dpi = 300, bg = "white")

# ═════════ PART 2: TCGA patient-level analysis ═════════
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
cdm <- as.data.frame(colData(se))[tum, ]; cdm$patient <- substr(colnames(se)[tum], 1, 12)

probes <- P$probe
bt <- assay(se[probes, tum]); colnames(bt) <- cdm$patient
meth <- colMeans(bt, na.rm = TRUE)                       # gene-level promoter beta (tumors)
bn <- assay(se[probes, nor]); meth_n <- colMeans(bn, na.rm = TRUE)   # normals

# probe-level tumor vs normal (includes cg23932332)
ptn <- data.frame(probe = probes,
                  mean_tumor = rowMeans(bt, na.rm = TRUE), mean_normal = rowMeans(bn, na.rm = TRUE)) %>%
  mutate(delta = mean_tumor - mean_normal)
write.csv(ptn, file.path(a, "19b_DUSP10_probe_tumor_vs_normal.csv"), row.names = FALSE)
cat("\n== DUSP10 probe-level tumor minus normal beta ==\n"); print(ptn, digits = 3)

vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens  <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row  <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
et   <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
en   <- which(substr(colnames(expr), 14, 15) == "11")
ex   <- expr[row, et]; names(ex) <- substr(colnames(expr)[et], 1, 12)
exn  <- unname(expr[row, en])

age_col <- intersect(c("age_at_index", "age_at_diagnosis"), colnames(cdm))[1]
age <- as.numeric(cdm[[age_col]]); if (median(age, na.rm = TRUE) > 200) age <- age / 365.25
dead <- tolower(cdm$vital_status) == "dead"
time <- as.numeric(ifelse(dead, cdm$days_to_death, cdm$days_to_last_follow_up)) / 30.44
st <- cdm$ajcc_pathologic_stage
stage_num <- ifelse(grepl("^Stage IV", st), 4, ifelse(grepl("^Stage III", st), 3,
                                                      ifelse(grepl("^Stage II", st), 2, ifelse(grepl("^Stage I", st), 1, NA))))
d <- data.frame(patient = cdm$patient, meth = meth[cdm$patient], expr = unname(ex[cdm$patient]),
                grade = factor(cdm$tumor_grade, c("G1", "G2", "G3", "G4")),
                time = time, event = as.integer(dead), age = age, stage_num = stage_num,
                stage = factor(ifelse(stage_num %in% 1:2, "Early", ifelse(stage_num %in% 3:4, "Late", NA)), c("Early", "Late")))
d$grade_num <- as.integer(d$grade)
d$grade_grp <- factor(ifelse(d$grade_num %in% 1:2, "Low", ifelse(d$grade_num %in% 3:4, "High", NA)), c("Low", "High"))
cat("\nTumors with expression:", sum(is.finite(d$expr)), "| with methylation:", sum(is.finite(d$meth)),
    "| normals (expr):", length(exn), "| normals (meth):", length(meth_n), "\n")

# helpers
two <- function(x, f) { k <- is.finite(x) & !is.na(f); l <- levels(f)
c(est = median(x[k & f == l[2]]) - median(x[k & f == l[1]]),
  p = wilcox.test(x[k & f == l[2]], x[k & f == l[1]])$p.value, n = sum(k)) }
sp <- function(x, y) { k <- is.finite(x) & is.finite(y)
ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE))
c(est = unname(ct$estimate), p = ct$p.value, n = sum(k)) }
tr <- function(feature, test, v) data.frame(feature, test, estimate = unname(v["est"]), p = unname(v["p"]), n = unname(v["n"]))

tests <- bind_rows(
  tr("Expression", "Tumor minus normal (mean VST)",
     c(est = mean(d$expr, na.rm = TRUE) - mean(exn, na.rm = TRUE), p = wilcox.test(d$expr, exn)$p.value, n = sum(is.finite(d$expr)) + length(exn))),
  tr("Expression", "Late minus Early stage (median VST)", two(d$expr, d$stage)),
  tr("Expression", "Stage I-IV ordered (Spearman rho)", sp(d$expr, d$stage_num)),
  tr("Expression", "High minus Low grade (median VST)", two(d$expr, d$grade_grp)),
  tr("Expression", "Grade G1-G4 ordered (Spearman rho)", sp(d$expr, d$grade_num)),
  tr("Gene methylation", "Tumor minus normal (delta beta)",
     c(est = mean(d$meth, na.rm = TRUE) - mean(meth_n, na.rm = TRUE), p = wilcox.test(d$meth, meth_n)$p.value, n = sum(is.finite(d$meth)) + length(meth_n))),
  tr("Gene methylation", "Late minus Early stage (median beta)", two(d$meth, d$stage)),
  tr("Gene methylation", "Stage I-IV ordered (Spearman rho)", sp(d$meth, d$stage_num)),
  tr("Gene methylation", "High minus Low grade (median beta)", two(d$meth, d$grade_grp)),
  tr("Gene methylation", "Grade G1-G4 ordered (Spearman rho)", sp(d$meth, d$grade_num)),
  tr("Methylation vs expression", "All tumors (Spearman rho)", sp(d$meth, d$expr)))
tests$p_BH <- p.adjust(tests$p, "BH")
write.csv(tests, file.path(a, "21_DUSP10_TCGA_tests.csv"), row.names = FALSE)
cat("\n== DUSP10 TCGA tests (BH over this table, exploratory) ==\n"); print(tests, digits = 3)

# existing DESeq2 result for comparison with script 38
cat("\n== Existing DESeq2 row (script 38) ==\n")
print(tryCatch(read.csv(file.path(d_dir, "nexus_stage_grade_dge.csv")) %>% filter(gene == G), error = function(e) conditionMessage(e)))
cat("\n== Existing survival row (script 45), compare with models below ==\n")
print(tryCatch(read.csv(file.path(a, "8_expr_vs_survival_gene.csv")) %>% filter(gene == G), error = function(e) conditionMessage(e)))

# survival models (HR per 1 SD)
d$expr_z <- as.numeric(scale(d$expr)); d$meth_z <- as.numeric(scale(d$meth))
okt <- is.finite(d$time) & d$time > 0
ds_e <- d[okt & is.finite(d$expr_z), ]
ds_b <- d[okt & is.finite(d$expr_z) & is.finite(d$meth_z), ]
cox_out <- function(f, dat, label) {
  fit <- coxph(f, data = dat); s <- summary(fit)
  data.frame(model = label, term = rownames(s$coefficients), HR = s$conf.int[, 1], lo = s$conf.int[, 3],
             hi = s$conf.int[, 4], p = s$coefficients[, 5], n = fit$n, events = fit$nevent, row.names = NULL) }
models <- rbind(
  cox_out(Surv(time, event) ~ expr_z, ds_e, "Expression alone"),
  cox_out(Surv(time, event) ~ expr_z + age + stage, ds_e, "Expression + age + stage"),
  cox_out(Surv(time, event) ~ meth_z, ds_b, "Methylation alone"),
  cox_out(Surv(time, event) ~ expr_z + meth_z, ds_b, "Expression + methylation"),
  cox_out(Surv(time, event) ~ expr_z + meth_z + age + stage, ds_b, "Expression + methylation + age + stage"))
write.csv(models, file.path(a, "20_DUSP10_survival_models.csv"), row.names = FALSE)
cat("\n== DUSP10 survival models (HR per 1 SD) ==\n"); print(models, digits = 3)

# figures: expression and methylation by group
mk <- function(val, nm) bind_rows(
  data.frame(panel = "Tumor vs normal", group = c(rep("Normal", sum(nm == "n")), rep("Tumor", sum(nm == "t"))), value = NA)[0, ])
fd <- bind_rows(
  data.frame(measure = "Expression (VST)", panel = "Tumor vs normal", group = "Normal", value = exn),
  data.frame(measure = "Expression (VST)", panel = "Tumor vs normal", group = "Tumor", value = d$expr),
  data.frame(measure = "Expression (VST)", panel = "Stage", group = ifelse(is.na(d$stage_num), NA, paste0("Stage ", as.roman(d$stage_num))), value = d$expr),
  data.frame(measure = "Expression (VST)", panel = "Grade", group = as.character(d$grade), value = d$expr),
  data.frame(measure = "Promoter methylation (beta)", panel = "Tumor vs normal", group = "Normal", value = meth_n),
  data.frame(measure = "Promoter methylation (beta)", panel = "Tumor vs normal", group = "Tumor", value = d$meth),
  data.frame(measure = "Promoter methylation (beta)", panel = "Stage", group = ifelse(is.na(d$stage_num), NA, paste0("Stage ", as.roman(d$stage_num))), value = d$meth),
  data.frame(measure = "Promoter methylation (beta)", panel = "Grade", group = as.character(d$grade), value = d$meth)) %>%
  filter(!is.na(group), is.finite(value))
fd$panel <- factor(fd$panel, c("Tumor vs normal", "Stage", "Grade"))
fB <- ggplot(fd, aes(group, value, fill = panel)) + geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  facet_grid(measure ~ panel, scales = "free", space = "free_x") +
  labs(title = "DUSP10 by tissue, stage and grade (TCGA-LIHC)", x = NULL, y = NULL) + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_DUSP10_expr_meth_by_group.png"), fB, width = 10, height = 6, dpi = 300, bg = "white")

ok <- is.finite(d$meth) & is.finite(d$expr)
fC <- ggplot(d[ok, ], aes(meth, expr, color = grade)) + geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE) +
  labs(title = "DUSP10: promoter methylation vs expression", x = "Mean promoter beta", y = "VST expression") + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "fig_DUSP10_meth_vs_expr.png"), fC, width = 6, height = 5, dpi = 300, bg = "white")

# ═════════ PART 3: GEO replication (helpers copied from script 50) ═════════
prep_annot <- function(f) { tab <- Table(getGEO(filename = f))
sc <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
list(ids = as.character(tab$ID), sy = lapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), trimws)) }
gene_vec <- function(ann, mat, gene) {
  idx <- which(vapply(ann$sy, function(x) gene %in% x, logical(1)))
  hit <- intersect(ann$ids[idx], rownames(mat)); if (length(hit) == 0) return(NULL)
  best <- hit[which.max(rowMeans(mat[hit, , drop = FALSE]))]
  list(x = mat[best, ], probe = best, nprobes = length(hit)) }
tnm_num <- function(x) { x <- toupper(trimws(as.character(x)))
ifelse(x == "I", 1, ifelse(x == "II", 2, ifelse(grepl("^III", x), 3, ifelse(grepl("^IV", x), 4, NA)))) }
test_gene <- function(d, mat, ann, gene, ds, linear) {
  gv <- gene_vec(ann, mat, gene)
  if (is.null(gv)) return(data.frame(dataset = ds, gene = gene, note = "no probe"))
  d$x <- unname(gv$x[d$id]); if (linear) d$x <- log2(pmax(d$x, 1))
  r <- data.frame(dataset = ds, gene = gene, probe = gv$probe, n_probes = gv$nprobes, note = "")
  tt <- d[!is.na(d$x) & !is.na(d$tumor), ]
  r$tn_diff <- mean(tt$x[tt$tumor]) - mean(tt$x[!tt$tumor])
  r$tn_p <- wilcox.test(tt$x[tt$tumor], tt$x[!tt$tumor])$p.value
  t <- tt[tt$tumor, ]
  t$stage_grp <- factor(ifelse(t$stage_num %in% 1:2, "Early", ifelse(t$stage_num %in% 3:4, "Late", NA)), c("Early", "Late"))
  s <- t[!is.na(t$time) & t$time > 0 & !is.na(t$event), ]; s$z <- as.numeric(scale(s$x))
  f1 <- tryCatch(summary(coxph(Surv(time, event) ~ z, data = s)), error = function(e) NULL)
  r$HR <- if (is.null(f1)) NA else f1$conf.int[1, 1]
  r$HR_low <- if (is.null(f1)) NA else f1$conf.int[1, 3]; r$HR_high <- if (is.null(f1)) NA else f1$conf.int[1, 4]
  r$cox_p <- if (is.null(f1)) NA else f1$coefficients[1, 5]
  f2 <- tryCatch(summary(coxph(Surv(time, event) ~ z + age + stage_grp, data = s)), error = function(e) NULL)
  r$HR_adj <- if (is.null(f2)) NA else f2$conf.int["z", 1]; r$cox_p_adj <- if (is.null(f2)) NA else f2$coefficients["z", 5]
  k <- t[!is.na(t$stage_num), ]
  ct <- suppressWarnings(cor.test(k$x, k$stage_num, method = "spearman", exact = FALSE))
  r$stage_rho <- unname(ct$estimate); r$stage_p <- ct$p.value
  e <- t[!is.na(t$stage_grp), ]
  r$late_early_diff <- median(e$x[e$stage_grp == "Late"]) - median(e$x[e$stage_grp == "Early"])
  r$late_early_p <- wilcox.test(e$x[e$stage_grp == "Late"], e$x[e$stage_grp == "Early"])$p.value
  r$n_surv <- nrow(s); r$events <- sum(s$event); r }

m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
tis <- tolower(trimws(cl$`Tissue Type`))
d1 <- data.frame(id = as.character(cl$Affy_GSM), tumor = ifelse(grepl("non", tis), FALSE, grepl("tumor", tis)),
                 time = as.numeric(cl$`Survival months`), event = as.numeric(cl$`Survival status`),
                 age = as.numeric(cl$Age), stage_num = tnm_num(cl$`TNM staging`))
d1 <- d1[d1$id %in% colnames(m1), ]
ann1 <- prep_annot(paste0(g, "GPL3921.annot.gz")); lin1 <- max(m1, na.rm = TRUE) > 100

m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]]); p2 <- p2[match(colnames(m2), rownames(p2)), ]
tis2 <- tolower(as.character(p2$`tissue:ch1`))
d2 <- data.frame(id = colnames(m2), tumor = ifelse(grepl("non", tis2), FALSE, grepl("tumor", tis2)),
                 time = as.numeric(as.character(p2$`duryears_os:ch1`)) * 12, event = as.numeric(as.character(p2$`event_os:ch1`)),
                 age = as.numeric(as.character(p2$`age (years):ch1`)), stage_num = tnm_num(p2$`tnm_staging_clinical:ch1`))
ann2 <- prep_annot(paste0(g, "GPL10558.annot.gz")); lin2 <- max(m2, na.rm = TRUE) > 100

r1 <- test_gene(d1, m1, ann1, G, "GSE14520", lin1)
r2 <- test_gene(d2, m2, ann2, G, "GSE76427", lin2)
geo <- bind_rows(r1, r2)
write.csv(geo, file.path(a, "22_DUSP10_GEO_tests.csv"), row.names = FALSE)
cat("\n== DUSP10 in GEO (same code as script 50) ==\n"); print(as.data.frame(geo), digits = 3)

# Kaplan-Meier: median split of expression (TCGA, GEO) and methylation (TCGA)
km_panel <- function(dd, nm) {
  dd <- dd[is.finite(dd$time) & dd$time > 0 & is.finite(dd$event) & is.finite(dd$x), ]
  dd$grp <- factor(ifelse(dd$x > median(dd$x), "High", "Low"), c("Low", "High"))
  lr <- survdiff(Surv(time, event) ~ grp, data = dd)
  plot(survfit(Surv(time, event) ~ grp, data = dd), col = c("#2C7BB6", "#D7301F"), lwd = 2,
       xlab = "Months", ylab = "Overall survival",
       main = sprintf("%s\nlog-rank p = %.3g (n = %d)", nm, 1 - pchisq(lr$chisq, 1), nrow(dd)), cex.main = 0.9)
  legend("bottomleft", c("Low", "High"), col = c("#2C7BB6", "#D7301F"), lwd = 2, bty = "n", cex = 0.8) }
k1 <- d1[d1$tumor, ]; gv1 <- gene_vec(ann1, m1, G); k1$x <- unname(gv1$x[k1$id]); if (lin1) k1$x <- log2(pmax(k1$x, 1))
k2 <- d2[d2$tumor, ]; gv2 <- gene_vec(ann2, m2, G); k2$x <- unname(gv2$x[k2$id]); if (lin2) k2$x <- log2(pmax(k2$x, 1))
png(file.path(fig_dir, "fig_DUSP10_km.png"), width = 2000, height = 2000, res = 250)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
km_panel(transform(d, x = expr), "DUSP10 expression, TCGA-LIHC")
km_panel(transform(d, x = meth), "DUSP10 promoter methylation, TCGA-LIHC")
km_panel(k1, "DUSP10 expression, GSE14520")
km_panel(k2, "DUSP10 expression, GSE76427")
dev.off()
cat("\nDone. Saved tables 19, 19b, 20, 21, 22 and 4 figures.\n")