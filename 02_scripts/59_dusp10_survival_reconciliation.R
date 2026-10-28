## 59: Reconcile the DUSP10 survival numbers
## Earlier report (Step 7):  HR 0.740, p = 0.0006, 366 patients, 132 deaths; stage-adjusted HR 0.789, p = 0.012
## Script 51 (table 23):     HR 0.77,  p = 0.003,  304 patients, 130 deaths
## Idea: script 51 took clinical data from the METHYLATION object and kept only time > 0.
## This script runs the same Cox model on (a) the expression object's clinical data and (b) the methylation
## object's clinical data, under three rules for building survival time, so we can see which choice gives which n.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(survival); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
G <- "DUSP10"

# ---- 1. DUSP10 expression, one primary tumor per patient ----
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ex <- setNames(as.numeric(expr[row, et]), substr(colnames(expr)[et], 1, 12))
cat("Primary tumors with expression:", length(ex), "\n")

# ---- 2. Clinical data from the two sources ----
need <- c("vital_status", "days_to_death", "days_to_last_follow_up")
cdv <- as.data.frame(colData(vsd))[et, ]; cdv$patient <- names(ex)
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cdm <- as.data.frame(colData(se))[tum, ]; cdm$patient <- substr(colnames(se)[tum], 1, 12)
cat("Clinical columns present in expression object:", all(need %in% names(cdv)),
    "| in methylation object:", all(need %in% names(cdm)), "\n")
stopifnot(all(need %in% names(cdv)), all(need %in% names(cdm)))

mk_df <- function(cd) {
  dead <- tolower(cd$vital_status) == "dead"
  age_col <- intersect(c("age_at_index", "age_at_diagnosis"), names(cd))[1]
  age <- if (is.na(age_col)) rep(NA_real_, nrow(cd)) else as.numeric(cd[[age_col]])
  if (median(age, na.rm = TRUE) > 200) age <- age / 365.25
  st <- if ("ajcc_pathologic_stage" %in% names(cd)) cd$ajcc_pathologic_stage else rep(NA_character_, nrow(cd))
  sn <- ifelse(grepl("^Stage IV", st), 4, ifelse(grepl("^Stage III", st), 3,
                                                 ifelse(grepl("^Stage II", st), 2, ifelse(grepl("^Stage I", st), 1, NA))))
  data.frame(patient = cd$patient, expr = unname(ex[cd$patient]), dead = dead,
             d_death = as.numeric(cd$days_to_death), d_fu = as.numeric(cd$days_to_last_follow_up), age = age,
             stage = factor(ifelse(sn %in% 1:2, "Early", ifelse(sn %in% 3:4, "Late", NA)), c("Early", "Late")))
}
src <- list(`expression object` = mk_df(cdv), `methylation object` = mk_df(cdm))

cat("\n== Why patients drop out: survival fields by source ==\n")
print(bind_rows(lapply(names(src), function(s) { d <- src[[s]]
data.frame(source = s, tumors = nrow(d), with_expr = sum(is.finite(d$expr)), deaths = sum(d$dead),
           deaths_missing_death_day = sum(d$dead & !is.finite(d$d_death)),
           alive_missing_followup = sum(!d$dead & !is.finite(d$d_fu)),
           followup_is_zero = sum(!d$dead & d$d_fu %in% 0 | d$dead & d$d_death %in% 0)) })))

# ---- 3. Cox model under different rules ----
cox_run <- function(d, source, rule) {
  time_days <- switch(rule,
                      "A: time>0 (script 51 rule)"         = ifelse(d$dead, d$d_death, d$d_fu),
                      "B: time>=0"                         = ifelse(d$dead, d$d_death, d$d_fu),
                      "C: dead w/o death day -> follow-up" = ifelse(d$dead & is.finite(d$d_death), d$d_death, d$d_fu))
  ok <- is.finite(time_days) & is.finite(d$expr) & (if (rule == "A: time>0 (script 51 rule)") time_days > 0 else time_days >= 0)
  s <- d[ok, ]; s$time <- time_days[ok] / 30.44; s$z <- as.numeric(scale(s$expr)); s$event <- as.integer(s$dead)
  f_sd   <- summary(coxph(Surv(time, event) ~ z, data = s))
  f_unit <- summary(coxph(Surv(time, event) ~ expr, data = s))
  f_st   <- tryCatch(summary(coxph(Surv(time, event) ~ z + stage, data = s)), error = function(e) NULL)
  f_as   <- tryCatch(summary(coxph(Surv(time, event) ~ z + age + stage, data = s)), error = function(e) NULL)
  data.frame(source = source, rule = rule, n = nrow(s), events = sum(s$event),
             HR_per_SD = f_sd$conf.int[1, 1], lo = f_sd$conf.int[1, 3], hi = f_sd$conf.int[1, 4], p = f_sd$coefficients[1, 5],
             HR_per_unit = f_unit$conf.int[1, 1], p_unit = f_unit$coefficients[1, 5],
             HR_stage_adj = if (is.null(f_st)) NA else f_st$conf.int["z", 1],
             p_stage_adj = if (is.null(f_st)) NA else f_st$coefficients["z", 5],
             HR_age_stage_adj = if (is.null(f_as)) NA else f_as$conf.int["z", 1],
             p_age_stage_adj = if (is.null(f_as)) NA else f_as$coefficients["z", 5], row.names = NULL) }

rules <- c("A: time>0 (script 51 rule)", "B: time>=0", "C: dead w/o death day -> follow-up")
res <- bind_rows(lapply(names(src), function(s) bind_rows(lapply(rules, function(r) cox_run(src[[s]], s, r)))))
write.csv(res, file.path(a, "29_DUSP10_survival_reconciliation.csv"), row.names = FALSE)
cat("\n== DUSP10 Cox model by clinical source and time rule ==\n"); print(as.data.frame(res), digits = 3)

cat("\n== Existing row from the earlier survival script (file 8) ==\n")
print(tryCatch(read.csv(file.path(a, "8_expr_vs_survival_gene.csv")) %>% filter(gene == G), error = function(e) conditionMessage(e)))

cat("\nHow to read it: find the row with n = 366 and 132 events. HR_per_unit or HR_per_SD there should match 0.740 (earlier report).\n")
cat("HR_stage_adj should match 0.789. The earlier report may use HR per 1 unit of VST, while table 23 uses HR per 1 SD.\n")
cat("Use ONE source and rule for DUSP10 everywhere (report, table 23) and state the HR scaling.\n")