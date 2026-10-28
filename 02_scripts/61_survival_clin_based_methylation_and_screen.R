## 61: Methylation survival (DUSP10 gene-level and cg23932332) and the 26-gene expression survival screen,
## all on the patient-level clinical table (same source and rule as scripts 16, 28 and 60; n = 366 for DUSP10).
## Reports HR per 1 SD AND per 1 unit, so the scaling is explicit.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(survival); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
key <- "cg23932332"

# ---- 1. Survival variables (as in script 16) ----
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

# generic Cox helper: returns HR, CI, p, n, events for the first term
run <- function(dat, xcol, extra = "") {
  f <- tryCatch(coxph(as.formula(paste("Surv(time, os_event) ~", xcol, extra)), data = dat), error = function(e) NULL)
  if (is.null(f)) return(c(HR = NA, lo = NA, hi = NA, p = NA, n = NA, events = NA))
  s <- summary(f); c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5], n = f$n, events = f$nevent) }

# ---- 2. 26-gene expression survival screen on the full patient set ----
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ids <- substr(colnames(expr)[et], 1, 12)
genes <- read.csv(file.path(a, "8_expr_vs_survival_gene.csv"))$gene
ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw <- match(ens, sub("\\..*$", "", rownames(expr)))
cat("Genes in screen:", length(genes), "| not found in expression matrix:", paste(genes[is.na(rw)], collapse = ", "), "\n")
genes <- genes[!is.na(rw)]; rw <- rw[!is.na(rw)]
keep <- ids %in% cs$submitter_id
sv <- cs[match(ids[keep], cs$submitter_id), ]
cat("Tumors with expression and survival:", nrow(sv), "| deaths:", sum(sv$os_event), "\n")

screen <- bind_rows(lapply(seq_along(genes), function(i) {
  d <- sv; d$x <- as.numeric(expr[rw[i], et][keep]); d$z <- as.numeric(scale(d$x))
  u <- run(d, "z"); su <- run(d, "x"); s <- run(d, "z", "+ stage"); as_ <- run(d, "z", "+ age + stage")
  data.frame(gene = genes[i], n = u["n"], events = u["events"], HR_per_SD = u["HR"], lo = u["lo"], hi = u["hi"], p = u["p"],
             HR_per_unit = su["HR"], HR_stage_adj = s["HR"], p_stage_adj = s["p"],
             HR_age_stage_adj = as_["HR"], p_age_stage_adj = as_["p"], row.names = NULL) }))
screen$padj_BH <- p.adjust(screen$p, "BH")
screen$padj_stage_adj_BH <- p.adjust(screen$p_stage_adj, "BH")
screen$padj_age_stage_adj_BH <- p.adjust(screen$p_age_stage_adj, "BH")
screen <- screen %>% arrange(p)
write.csv(screen, file.path(a, "31_expr_survival_screen_clin_based.csv"), row.names = FALSE)
cat("\n== 26-gene expression survival screen (full patient set) ==\n")
print(screen %>% dplyr::select(gene, n, events, HR_per_SD, HR_per_unit, p, padj_BH, HR_stage_adj, p_stage_adj, padj_stage_adj_BH), digits = 3)
cat("\nEarlier report: DUSP10 padj 0.016, NPC1 0.044, HDAC11 0.060 (those came from the older run).\n")

# ---- 3. Methylation survival: DUSP10 gene-level and cg23932332 ----
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
probes <- read.csv(file.path(a, "19_DUSP10_probe_table.csv"))$probe
M <- assay(se[probes, tum])
md <- data.frame(submitter_id = substr(colnames(se)[tum], 1, 12),
                 gene_meth = colMeans(M, na.rm = TRUE), key_meth = as.numeric(assay(se[key, tum])))
d <- inner_join(md, cs, by = "submitter_id")
ed <- data.frame(submitter_id = ids, dusp10 = as.numeric(expr[rw[match("DUSP10", genes)], et]))
d <- left_join(d, ed, by = "submitter_id")
d$z_gene <- as.numeric(scale(d$gene_meth)); d$z_key <- as.numeric(scale(d$key_meth)); d$z_expr <- as.numeric(scale(d$dusp10))
cat("\nMethylation tumors matched to survival:", nrow(d), "| deaths:", sum(d$os_event), "\n")

mrow <- function(label, v) data.frame(model = label, t(v), row.names = NULL)
meth <- bind_rows(
  mrow("Gene-level promoter methylation, HR per 1 SD",                run(d, "z_gene")),
  mrow("Gene-level promoter methylation, HR per 1 unit beta",         run(d, "gene_meth")),
  mrow("Gene-level promoter methylation + age + stage, per 1 SD",     run(d, "z_gene", "+ age + stage")),
  mrow("cg23932332, HR per 1 SD",                                     run(d, "z_key")),
  mrow("cg23932332, HR per 1 unit beta",                              run(d, "key_meth")),
  mrow("cg23932332 + age + stage, per 1 SD",                          run(d, "z_key", "+ age + stage")),
  mrow("cg23932332 + DUSP10 expression, per 1 SD (methylation term)", run(d, "z_key", "+ z_expr")),
  mrow("DUSP10 expression + cg23932332, per 1 SD (expression term)",  run(d, "z_expr", "+ z_key")))
write.csv(meth, file.path(a, "31b_methylation_survival_clin_based.csv"), row.names = FALSE)
cat("\n== DUSP10 methylation and survival (full patient set) ==\n"); print(meth, digits = 3)
cat("\nThe probe-level padj in file 10 was adjusted over many probes, so only the raw p-values here are comparable to table 23.\n")
cat("Earlier table values on the 304-patient set: gene-level HR 1.03 (p = 0.73); cg23932332 HR 1.09 (padj 0.71, scaling unknown).\n")
cat("To see how script 46 scaled the probe HR: grep -n 'scale\\|coxph' 02_scripts/46_methylation_survival.R\n")