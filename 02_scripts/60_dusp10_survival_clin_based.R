## 60: DUSP10 survival using the patient-level clinical table (the same source and rule as scripts 16 and 28)
## Goal: reproduce the earlier report (n = 366, 132 deaths, HR 0.740, p = 0.0006; stage-adjusted HR 0.789, p = 0.012).
## Scripts 45/46/48/51 took clinical data from the methylation/expression objects, which lack follow-up times for ~60 living patients (n = 304).
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(survival); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
G <- "DUSP10"

# ---- 1. Survival variables, exactly as in script 16 ----
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")
cat("Clinical table:", nrow(clin), "patients x", ncol(clin), "columns\n")
clin$os_event <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$os_time  <- ifelse(!is.na(clin$days_to_death), clin$days_to_death, clin$days_to_last_follow_up)
cs <- clin[!is.na(clin$os_time) & clin$os_time > 0, ]
cat("Patients with survival time > 0:", nrow(cs), "| deaths:", sum(cs$os_event), "\n")

stage_src <- intersect(c("ajcc_pathologic_stage", "tumor_stage"), names(clin))[1]
st <- as.character(clin[[stage_src]])
cat("Stage column used:", stage_src, "| values:", paste(head(unique(st), 8), collapse = " | "), "\n")
sn <- ifelse(grepl("(?i)stage iv", st, perl = TRUE), 4, ifelse(grepl("(?i)stage iii", st, perl = TRUE), 3,
                                                               ifelse(grepl("(?i)stage ii", st, perl = TRUE), 2, ifelse(grepl("(?i)stage i", st, perl = TRUE), 1, NA))))
clin$stage_num <- sn
clin$stage <- factor(ifelse(sn %in% 1:2, "Early", ifelse(sn %in% 3:4, "Late", NA)), c("Early", "Late"))
if ("age_at_index" %in% names(clin)) clin$age <- as.numeric(clin$age_at_index) else
  clin$age <- abs(as.numeric(clin$days_to_birth)) / 365.25
cs <- clin[!is.na(clin$os_time) & clin$os_time > 0, ]

# ---- 2. DUSP10 expression, one primary tumor per patient ----
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ex <- data.frame(submitter_id = substr(colnames(expr)[et], 1, 12), expr = as.numeric(expr[row, et]))

d <- inner_join(ex, cs[, c("submitter_id", "os_time", "os_event", "stage", "stage_num", "age")], by = "submitter_id")
d$time <- d$os_time / 30.44; d$z <- as.numeric(scale(d$expr))
cat("\nTumors with expression:", nrow(ex), "| matched to survival:", nrow(d), "| deaths:", sum(d$os_event), "\n")

# ---- 3. Cox models ----
get <- function(f, term) { s <- summary(f); c(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3], hi = s$conf.int[term, 4], p = s$coefficients[term, 5], n = f$n, events = f$nevent) }
m <- list(
  `Expression alone, HR per 1 SD`          = get(coxph(Surv(time, os_event) ~ z, data = d), "z"),
  `Expression alone, HR per 1 VST unit`    = get(coxph(Surv(time, os_event) ~ expr, data = d), "expr"),
  `+ stage (early/late), HR per 1 SD`      = get(coxph(Surv(time, os_event) ~ z + stage, data = d), "z"),
  `+ stage (I-IV numeric), HR per 1 SD`    = get(coxph(Surv(time, os_event) ~ z + stage_num, data = d), "z"),
  `+ age + stage (early/late), per 1 SD`   = get(coxph(Surv(time, os_event) ~ z + age + stage, data = d), "z"))
res <- bind_rows(lapply(names(m), function(k) data.frame(model = k, t(m[[k]]), row.names = NULL)))
write.csv(res, file.path(a, "30_DUSP10_survival_clin_based.csv"), row.names = FALSE)
cat("\n== DUSP10 survival, patient-level clinical table ==\n"); print(res, digits = 3)

cat("\nTargets from the earlier report: n = 366, deaths = 132, HR 0.740 (p = 0.0006), stage-adjusted HR 0.789 (p = 0.012).\n")
cat("If n and deaths match, the 366 definition is reproduced. Then check whether 0.740 appears per 1 SD or per 1 VST unit.\n")
cat("Script 16 filters os_time > 0 and uses days_to_death when it is not missing, else days_to_last_follow_up.\n")