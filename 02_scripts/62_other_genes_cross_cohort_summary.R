## 62: Cross-cohort summary for other genes (SMAD2, MYC, MAPK8, HDAC11), same layout as file 23 (DUSP10)
## One row per question, one column per cohort: TCGA-LIHC, GSE14520, GSE76427 (expression), GSE56588, GSE54503 (methylation).
## Everything for TCGA is recomputed here from the raw objects, using the SAME rules as the DUSP10 work:
##   - survival from the patient-level clinical table (script 16/60/61 rule), HR per 1 SD
##   - methylation = mean beta over the gene's probes in meth_expr_corr_probe_level.csv
##   - "strongest probe" = the gene's probe with the largest |rho| with expression in that file
## Raw p-values everywhere (no padj). GSE56588 and GSE54503 have methylation only; GEO expression cohorts have no methylation.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(GEOquery); library(survival); library(dplyr); library(tidyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; g <- "data/raw/geo/"
genes <- c("SMAD2", "MYC", "MAPK8", "HDAC11")   # add other symbols here (needs probes in meth_expr_corr_probe_level.csv)

mk <- function(layer, analysis, dataset, estimate, p, lo = NA, hi = NA, stat, n = NA) {
  data.frame(layer = layer, analysis = analysis, dataset = dataset, statistic = stat, estimate = estimate,
             ci_low = lo, ci_high = hi, p = p, n = n, stringsAsFactors = FALSE) }

# ================= 1. Survival base: patient-level clinical table =================
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")
clin$os_event <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$os_time  <- ifelse(!is.na(clin$days_to_death), clin$days_to_death, clin$days_to_last_follow_up)
st <- as.character(clin$ajcc_pathologic_stage)
sn <- ifelse(grepl("(?i)stage iv", st, perl = TRUE), 4, ifelse(grepl("(?i)stage iii", st, perl = TRUE), 3,
                                                               ifelse(grepl("(?i)stage ii", st, perl = TRUE), 2, ifelse(grepl("(?i)stage i", st, perl = TRUE), 1, NA))))
clin$stage <- factor(ifelse(sn %in% 1:2, "Early", ifelse(sn %in% 3:4, "Late", NA)), c("Early", "Late"))
clin$age <- if ("age_at_index" %in% names(clin)) as.numeric(clin$age_at_index) else abs(as.numeric(clin$days_to_birth)) / 365.25
cs <- clin[!is.na(clin$os_time) & clin$os_time > 0, c("submitter_id", "os_time", "os_event", "stage", "age")]
cs$time <- cs$os_time / 30.44
run <- function(dat, xcol, extra = "") {
  f <- tryCatch(coxph(as.formula(paste("Surv(time, os_event) ~", xcol, extra)), data = dat), error = function(e) NULL)
  if (is.null(f)) return(c(HR = NA, lo = NA, hi = NA, p = NA, n = NA))
  s <- summary(f); c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5], n = f$n) }

# ================= 2. TCGA expression and methylation objects =================
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
bcE <- substr(colnames(expr), 14, 15)
et <- which(bcE == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
en <- which(bcE == "11"); ids_e <- substr(colnames(expr)[et], 1, 12)

se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
cdm <- as.data.frame(colData(se))[tum, ]; cdm$patient <- substr(colnames(se)[tum], 1, 12)
stx <- cdm$ajcc_pathologic_stage
stn <- ifelse(grepl("^Stage IV", stx), 4, ifelse(grepl("^Stage III", stx), 3,
                                                 ifelse(grepl("^Stage II", stx), 2, ifelse(grepl("^Stage I", stx), 1, NA))))
cdm$stage_num <- stn
cdm$stage <- factor(ifelse(stn %in% 1:2, "Early", ifelse(stn %in% 3:4, "Late", NA)), c("Early", "Late"))
gr <- factor(cdm$tumor_grade, c("G1", "G2", "G3", "G4")); cdm$grade_num <- as.integer(gr)
cdm$grade_grp <- factor(ifelse(cdm$grade_num %in% 1:2, "Low", ifelse(cdm$grade_num %in% 3:4, "High", NA)), c("Low", "High"))

two <- function(x, f) { k <- is.finite(x) & !is.na(f); l <- levels(f)
c(est = median(x[k & f == l[2]]) - median(x[k & f == l[1]]),
  p = wilcox.test(x[k & f == l[2]], x[k & f == l[1]])$p.value, n = sum(k)) }
sp <- function(x, y) { k <- is.finite(x) & is.finite(y)
ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE))
c(est = unname(ct$estimate), p = ct$p.value, n = sum(k)) }

pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))

tcga_rows <- function(G, probes, key) {
  tc <- "TCGA-LIHC"
  ens <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
  row <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
  ex_t <- setNames(as.numeric(expr[row, et]), ids_e); ex_n <- as.numeric(expr[row, en])
  gm_t <- setNames(colMeans(assay(se[probes, tum]), na.rm = TRUE), cdm$patient)
  gm_n <- colMeans(assay(se[probes, nor]), na.rm = TRUE)
  kb_t <- setNames(as.numeric(assay(se[key, tum])), cdm$patient); kb_n <- as.numeric(assay(se[key, nor]))
  d <- cdm[, c("patient", "stage", "stage_num", "grade_num", "grade_grp")]
  d$expr <- unname(ex_t[d$patient]); d$meth <- unname(gm_t[d$patient]); d$kb <- unname(kb_t[d$patient])
  r <- function(layer, an, stat, v, n = NA, lo = NA, hi = NA) mk(layer, an, tc, unname(v["est"]), unname(v["p"]), lo, hi, stat, unname(v["n"]))
  wt <- function(x, y) c(est = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE), p = wilcox.test(x, y)$p.value, n = sum(is.finite(x)) + sum(is.finite(y)))
  out <- bind_rows(
    r("Expression", "Tumor vs normal", "mean difference", wt(d$expr, ex_n)),
    r("Expression", "Late vs early stage", "median difference", two(d$expr, d$stage)),
    r("Expression", "Stage I-IV trend", "Spearman rho", sp(d$expr, d$stage_num)),
    r("Expression", "High vs low grade", "median difference", two(d$expr, d$grade_grp)),
    r("Expression", "Grade G1-G4 trend", "Spearman rho", sp(d$expr, d$grade_num)),
    r("Methylation", "Tumor vs normal", "delta beta", wt(d$meth, gm_n)),
    r("Methylation", "Late vs early stage", "median difference", two(d$meth, d$stage)),
    r("Methylation", "Grade G1-G4 trend", "Spearman rho", sp(d$meth, d$grade_num)),
    r("Methylation", "Methylation vs expression", "Spearman rho", sp(d$meth, d$expr)))
  pk <- paste0("Methylation, probe ", key, " (strongest rho)")
  out <- bind_rows(out,
                   r(pk, "Tumor vs normal", "delta beta", wt(d$kb, kb_n)),
                   r(pk, "Methylation vs expression", "Spearman rho", sp(d$kb, d$expr)),
                   r(pk, "Grade G1-G4 trend", "Spearman rho", sp(d$kb, d$grade_num)),
                   r(pk, "Stage I-IV trend", "Spearman rho", sp(d$kb, d$stage_num)))
  # survival (patient-level clinical table), HR per 1 SD
  sv <- inner_join(data.frame(submitter_id = ids_e, x = unname(ex_t)), cs, by = "submitter_id")
  sv$z <- as.numeric(scale(sv$x))
  u <- run(sv, "z"); ad <- run(sv, "z", "+ age + stage")
  mv <- inner_join(data.frame(submitter_id = cdm$patient, m = unname(gm_t[cdm$patient]), k = unname(kb_t[cdm$patient])), cs, by = "submitter_id")
  mv$zm <- as.numeric(scale(mv$m)); mv$zk <- as.numeric(scale(mv$k))
  um <- run(mv, "zm"); uk <- run(mv, "zk")
  bind_rows(out,
            mk("Expression", "Survival, unadjusted", tc, u["HR"], u["p"], u["lo"], u["hi"], "HR per 1 SD", u["n"]),
            mk("Expression", "Survival, adjusted (age + stage)", tc, ad["HR"], ad["p"], ad["lo"], ad["hi"], "HR per 1 SD", ad["n"]),
            mk("Methylation", "Survival, unadjusted", tc, um["HR"], um["p"], um["lo"], um["hi"], "HR per 1 SD", um["n"]),
            mk(pk, "Survival, unadjusted", tc, uk["HR"], uk["p"], uk["lo"], uk["hi"], "HR per 1 SD", uk["n"]))
}

# ================= 3. GEO expression cohorts (helpers copied from script 51) =================
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
  if (is.null(gv)) return(NULL)
  d$x <- unname(gv$x[d$id]); if (linear) d$x <- log2(pmax(d$x, 1))
  r <- data.frame(dataset = ds, gene = gene)
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
  r$n_surv <- nrow(s); r }

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

geo_rows <- function(r) { if (is.null(r)) return(NULL); bind_rows(
  mk("Expression", "Tumor vs normal", r$dataset, r$tn_diff, r$tn_p, stat = "mean difference"),
  mk("Expression", "Late vs early stage", r$dataset, r$late_early_diff, r$late_early_p, stat = "median difference"),
  mk("Expression", "Stage I-IV trend", r$dataset, r$stage_rho, r$stage_p, stat = "Spearman rho"),
  mk("Expression", "Survival, unadjusted", r$dataset, r$HR, r$cox_p, r$HR_low, r$HR_high, stat = "HR per 1 SD", n = r$n_surv),
  mk("Expression", "Survival, adjusted (age + stage)", r$dataset, r$HR_adj, r$cox_p_adj, stat = "HR per 1 SD", n = r$n_surv)) }

# ================= 4. GSE56588 (from files 17 / 17b) and GSE54503 (paired, cached series matrix) =================
f17 <- file.path(a, "17_GSE56588_meth_7genes.csv"); f17b <- file.path(a, "17b_GSE56588_probe_delta.csv")
t17 <- if (file.exists(f17)) read.csv(f17) else NULL
t17b <- if (file.exists(f17b)) read.csv(f17b) else NULL

B54 <- NULL
f54 <- paste0(g, "GSE54503_series_matrix.txt.gz")
if (file.exists(f54)) {
  e54 <- getGEO(filename = f54, getGPL = FALSE); p54 <- pData(e54); B <- exprs(e54)
  chr <- grep(":ch1$", names(p54), value = TRUE)
  txt <- tolower(apply(p54[, c("title", chr), drop = FALSE], 1, paste, collapse = " "))
  p54$tissue <- ifelse(grepl("non-?tumou?r|adjacent|normal|benign|surrounding", txt), "Normal",
                       ifelse(grepl("tumou?r|hcc|carcinoma", txt), "Tumor", NA))
  p54 <- p54[order(rownames(p54)), ]
  ni <- which(p54$tissue == "Normal"); ok <- ni < nrow(p54) & p54$tissue[pmin(ni + 1, nrow(p54))] == "Tumor"
  nor_g <- rownames(p54)[ni[ok]]; tum_g <- rownames(p54)[ni[ok] + 1]
  rs <- grep("^rs", rownames(B), value = TRUE)
  same <- sapply(seq_along(nor_g), function(i) cor(B[rs, nor_g[i]], B[rs, tum_g[i]], use = "complete.obs"))
  other <- sapply(seq_along(nor_g), function(i) cor(B[rs, nor_g[i]], B[rs, tum_g[(i %% length(tum_g)) + 1]], use = "complete.obs"))
  if (median(same) > median(other) + 0.2) { B54 <- B; cat("GSE54503: pairing confirmed by genotype probes;", length(nor_g), "pairs.\n")
  } else cat("GSE54503 skipped: genotype check did not confirm the pairing.\n")
} else cat("GSE54503 series matrix not found in", g, ": GSE54503 column will be empty.\n")
paired <- function(x, y) c(est = mean(x - y, na.rm = TRUE), p = wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value)

# ================= 5. Loop over genes =================
fmt <- function(est, ci_low, ci_high, p, analysis) {
  plab <- ifelse(grepl("padj", analysis), "padj", "p")
  ifelse(is.na(ci_low),
         ifelse(is.na(p), sprintf("%.3g", est), sprintf("%.3g (%s=%.2g)", est, plab, p)),
         sprintf("%.2f [%.2f-%.2f] (%s=%.2g)", est, ci_low, ci_high, plab, p)) }
lev <- c("TCGA-LIHC", "GSE14520", "GSE76427", "GSE56588", "GSE54503")
all_long <- list()
for (G in genes) {
  cat("\n==================== ", G, " ====================\n", sep = "")
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No probes for", G, "in meth_expr_corr_probe_level.csv: skipped.\n"); next }
  probes <- pg$probe; key <- pg$probe[which.max(abs(pg$rho))]
  cat("Probes:", length(probes), "| strongest rho with expression:", key, "(rho =", round(pg$rho[pg$probe == key], 3), ")\n")
  pk <- paste0("Methylation, probe ", key, " (strongest rho)")
  rows <- list(tcga_rows(G, probes, key),
               geo_rows(test_gene(d1, m1, ann1, G, "GSE14520", lin1)),
               geo_rows(test_gene(d2, m2, ann2, G, "GSE76427", lin2)))
  if (!is.null(t17)) { x <- t17[t17$gene == G, ]
  if (nrow(x) == 1) rows <- c(rows, list(mk("Methylation", "Tumor vs normal", "GSE56588", x$delta_GSE56588, x$p_wilcox, stat = "delta beta"))) }
  if (!is.null(t17b)) { x <- t17b[t17b$gene == G & t17b$probe == key, ]
  if (nrow(x) == 1) rows <- c(rows, list(mk(pk, "Tumor vs normal", "GSE56588", x$delta, NA, stat = "delta beta"))) }
  if (!is.null(B54)) {
    pb <- intersect(probes, rownames(B54))
    if (length(pb) > 0) { gt <- colMeans(B54[pb, tum_g, drop = FALSE], na.rm = TRUE); gn <- colMeans(B54[pb, nor_g, drop = FALSE], na.rm = TRUE)
    r1 <- paired(gt, gn); rows <- c(rows, list(mk("Methylation", "Tumor vs normal", "GSE54503", r1["est"], r1["p"], stat = "delta beta", n = length(tum_g)))) }
    if (key %in% rownames(B54)) { r2 <- paired(B54[key, tum_g], B54[key, nor_g])
    rows <- c(rows, list(mk(pk, "Tumor vs normal", "GSE54503", r2["est"], r2["p"], stat = "delta beta", n = length(tum_g)))) }
  }
  long <- bind_rows(rows) %>% mutate(gene = G, dataset = factor(dataset, lev)) %>% arrange(layer, analysis, dataset)
  wide <- long %>% mutate(txt = fmt(estimate, ci_low, ci_high, p, analysis)) %>%
    dplyr::select(layer, analysis, statistic, dataset, txt) %>% pivot_wider(names_from = dataset, values_from = txt)
  write.csv(wide, file.path(a, paste0("32_", G, "_cross_cohort_summary.csv")), row.names = FALSE)
  print(as.data.frame(wide), right = FALSE)
  all_long[[G]] <- long
}
write.csv(bind_rows(all_long), file.path(a, "32_other_genes_cross_cohort_long.csv"), row.names = FALSE)
cat("\nSaved 32_<GENE>_cross_cohort_summary.csv for each gene and 32_other_genes_cross_cohort_long.csv.\n")
cat("Notes: p-values are raw (no padj). TCGA expression is VST, GEO is log2: compare direction, not size.\n")
cat("Survival: HR per 1 SD, TCGA from the patient-level clinical table. NA = analysis not possible in that cohort.\n")
cat("The strongest probe is chosen by |rho| with expression in TCGA, so its TCGA association values are optimistic (selected on that same data).\n")