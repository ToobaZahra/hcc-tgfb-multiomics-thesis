## 45: expression vs overall survival (Cox per gene + median-split Kaplan-Meier), 26 nexus genes, tumors only
## separate analysis (expression vs grade = table 1, vs stage = table 2, vs survival = this script)
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(survival)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; fig_dir <- "04_figures"

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Clinical data: overall survival (OS), age, stage ──
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
cdm <- as.data.frame(colData(se))
cdm$patient <- substr(colnames(se), 1, 12)
cdm <- cdm[substr(colnames(se), 14, 15) == "01", ]; cdm <- cdm[!duplicated(cdm$patient), ]
cat("Candidate survival/age columns:\n"); print(grep("vital|days|age|death|follow|survival", colnames(cdm), value = TRUE, ignore.case = TRUE))
need <- c("vital_status", "days_to_death", "days_to_last_follow_up")
if (length(setdiff(need, colnames(cdm))) > 0)
  stop("Missing columns: ", paste(setdiff(need, colnames(cdm)), collapse = ", "),
       " -- use tcga_lihc_clinical.rds instead (submitter_id key)")

age_col <- intersect(c("age_at_index", "age_at_diagnosis"), colnames(cdm))[1]
age <- if (is.na(age_col)) NA_real_ else as.numeric(cdm[[age_col]])
if (!is.na(age_col) && median(age, na.rm = TRUE) > 200) age <- age / 365.25      # days -> years
dead <- tolower(cdm$vital_status) == "dead"
time <- as.numeric(ifelse(dead, cdm$days_to_death, cdm$days_to_last_follow_up)) / 30.44   # months
stage <- ifelse(grepl("^Stage I{1,2}[ABC]?$", cdm$ajcc_pathologic_stage), "Early",
                ifelse(grepl("^Stage (III|IV)", cdm$ajcc_pathologic_stage), "Late", NA))
clin <- data.frame(patient = cdm$patient, time = time, event = as.integer(dead),
                   age = age, stage = factor(stage, c("Early", "Late")))
clin <- clin[!is.na(clin$time) & clin$time > 0, ]
cat("Patients with survival data:", nrow(clin), "| deaths:", sum(clin$event), "\n")

# ── 2. VST expression (tumors) ───────────────────────────
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
idx <- match(sym2ens, sub("\\..*$", "", rownames(expr))); ok <- !is.na(idx)
ex <- expr[idx[ok], , drop = FALSE]; rownames(ex) <- names(sym2ens)[ok]
et <- which(substr(colnames(ex), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(ex)[et], 1, 12))]
ex <- ex[, et, drop = FALSE]; colnames(ex) <- substr(colnames(ex), 1, 12)
clin <- clin[clin$patient %in% colnames(ex), ]
ex <- ex[, clin$patient, drop = FALSE]
cat("Patients used:", nrow(clin), "| deaths:", sum(clin$event), "\n")

# ── 3. Cox per gene (HR per 1 SD) + median-split log-rank ──
cox_hr <- function(fit) { s <- summary(fit); c(HR = s$conf.int["x", 1], lo = s$conf.int["x", 3],
                                               hi = s$conf.int["x", 4], p = s$coefficients["x", 5]) }
res <- do.call(rbind, lapply(rownames(ex), function(g) {
  d <- clin; d$x <- as.numeric(scale(ex[g, ]))
  u <- cox_hr(coxph(Surv(time, event) ~ x, data = d))
  ad <- tryCatch(cox_hr(coxph(Surv(time, event) ~ x + age + stage, data = d)),
                 error = function(e) c(HR = NA, lo = NA, hi = NA, p = NA))
  grp <- ifelse(ex[g, ] > median(ex[g, ]), "High", "Low")
  lr <- survdiff(Surv(time, event) ~ grp, data = d)
  data.frame(gene = g, n = nrow(d), events = sum(d$event),
             HR = u["HR"], HR_low = u["lo"], HR_high = u["hi"], p_cox = u["p"],
             HR_adj = ad["HR"], p_adj_model = ad["p"], logrank_p = 1 - pchisq(lr$chisq, 1), row.names = NULL)
}))
res$padj_cox <- p.adjust(res$p_cox, "BH")
res$padj_adj_model <- p.adjust(res$p_adj_model, "BH")
res$padj_logrank <- p.adjust(res$logrank_p, "BH")
res$survival_direction <- ifelse(res$padj_cox < 0.05,
                                 ifelse(res$HR > 1, "High expression = worse survival", "High expression = better survival"),
                                 "No association")

# ── 4. Save (survival only; grade and stage are in tables 1 and 2) ──
res <- res %>% arrange(p_cox)
write.csv(res, file.path(a, "8_expr_vs_survival_gene.csv"), row.names = FALSE)
print(as.data.frame(res %>% dplyr::select(gene, HR, HR_low, HR_high, p_cox, padj_cox, HR_adj, padj_adj_model,
                                          padj_logrank, survival_direction)), digits = 3)

# ── 5. Kaplan-Meier for the 6 strongest genes ────────────
top <- head(res$gene, 6)
png(file.path(fig_dir, "fig_km_top_genes.png"), width = 2400, height = 1600, res = 250)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (g in top) {
  d <- clin; d$grp <- factor(ifelse(ex[g, ] > median(ex[g, ]), "High", "Low"), c("Low", "High"))
  plot(survfit(Surv(time, event) ~ grp, data = d), col = c("#2C7BB6", "#D7301F"), lwd = 2,
       xlab = "Months", ylab = "Overall survival",
       main = sprintf("%s (log-rank p = %.3g)", g, res$logrank_p[res$gene == g]))
  legend("bottomleft", c("Low", "High"), col = c("#2C7BB6", "#D7301F"), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()