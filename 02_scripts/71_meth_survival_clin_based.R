## 71: Methylation vs overall survival on the patient-level clinical table (same source and rule as scripts 16, 28, 60, 61, 62)
## Replaces the clinical source of script 46 (methylation colData, which lacks follow-up for ~60 living patients).
## Everything else is copied from script 46: same 26 genes, same promoter definition (TSS -1500/+500), same Cox models
## (HR per 1 SD, + age + stage), same median-split log-rank, same BH adjustment, SAME OUTPUT COLUMNS as files 9 and 10,
## so script 47 only needs its two file names changed.
## Outputs (03_results/association):
##   41_meth_vs_survival_gene_clin_based.csv    (replaces file 9)
##   42_meth_vs_survival_probe_clin_based.csv   (replaces file 10)
## Figure: 04_figures/fig_km_methylation_top_clin_based.png
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(survival)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; fig_dir <- "04_figures"

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Promoter probes of the 26 genes (identical to script 46) ──
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit  <- vapply(g_list, function(g) any(g %in% genes), logical(1))
keep <- hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)
cand <- ann[keep, ]; cand$gene_list <- g_list[keep]
ent <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")); ent <- ent[!is.na(ent)]
tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")[ent]
prom <- unlist(promoters(tx_by_gene, upstream = 1500, downstream = 500))
prom$symbol <- names(ent)[match(names(prom), ent)]
prom <- keepStandardChromosomes(prom, pruning.mode = "coarse")
probe_gr <- GRanges(cand$chrm_A, IRanges(cand$beg_A, cand$beg_A + 1))
ov <- suppressWarnings(findOverlaps(probe_gr, prom))
probe_gene <- unique(data.frame(probe = cand$probe[queryHits(ov)], gene = prom$symbol[subjectHits(ov)]))
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]

# ── 2. Survival variables from the patient-level clinical table (same rule as scripts 16 and 61) ──
clin0 <- readRDS("data/raw/tcga_lihc_clinical.rds")
clin0$os_event <- ifelse(clin0$vital_status == "Dead", 1, 0)
clin0$os_time  <- ifelse(!is.na(clin0$days_to_death), clin0$days_to_death, clin0$days_to_last_follow_up)
st <- as.character(clin0$ajcc_pathologic_stage)
sn <- ifelse(grepl("(?i)stage iv", st, perl = TRUE), 4, ifelse(grepl("(?i)stage iii", st, perl = TRUE), 3,
                                                               ifelse(grepl("(?i)stage ii", st, perl = TRUE), 2, ifelse(grepl("(?i)stage i", st, perl = TRUE), 1, NA))))
stage <- factor(ifelse(sn %in% 1:2, "Early", ifelse(sn %in% 3:4, "Late", NA)), c("Early", "Late"))
age <- if ("age_at_index" %in% names(clin0)) as.numeric(clin0$age_at_index) else abs(as.numeric(clin0$days_to_birth)) / 365.25
clin <- data.frame(patient = clin0$submitter_id, time = clin0$os_time / 30.44, event = clin0$os_event,
                   age = age, stage = stage)
clin <- clin[!is.na(clin$time) & clin$time > 0, ]
clin <- clin[!duplicated(clin$patient), ]
cat("Patients with survival data (clinical table):", nrow(clin), "| deaths:", sum(clin$event), "\n")

# ── 3. Beta matrices: primary tumors, one sample per patient, matched to the clinical table ──
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
bt <- assay(se[unique(probe_gene$probe), tum]); colnames(bt) <- substr(colnames(se)[tum], 1, 12)
clin <- clin[clin$patient %in% colnames(bt), ]
bt <- bt[, clin$patient, drop = FALSE]
cat("Tumors with methylation AND survival:", nrow(clin), "| deaths:", sum(clin$event), "\n")
gb <- t(sapply(unique(probe_gene$gene), function(gn)
  colMeans(bt[probe_gene$probe[probe_gene$gene == gn], , drop = FALSE], na.rm = TRUE)))

# ── 4. Cox helpers (identical to script 46) ──
cox_x <- function(x, adj = FALSE) {
  d <- clin; d$x <- as.numeric(scale(x)); d <- d[is.finite(d$x), ]
  if (nrow(d) < 50 || sum(d$event) < 10) return(NULL)
  f <- if (adj) Surv(time, event) ~ x + age + stage else Surv(time, event) ~ x
  fit <- tryCatch(coxph(f, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  c(n = fit$n, events = fit$nevent, HR = s$conf.int["x", 1], lo = s$conf.int["x", 3],
    hi = s$conf.int["x", 4], p = s$coefficients["x", 5])
}
cox_table <- function(M) {
  out <- do.call(rbind, lapply(rownames(M), function(f) {
    x <- M[f, ]; u <- cox_x(x); if (is.null(u)) return(NULL)
    ad <- cox_x(x, TRUE)
    k <- is.finite(x); d <- clin[k, ]; d$grp <- ifelse(x[k] > median(x[k]), "High", "Low")
    lr <- tryCatch(survdiff(Surv(time, event) ~ grp, data = d), error = function(e) NULL)
    data.frame(feature = f, mean_beta = mean(x, na.rm = TRUE), n = u["n"], events = u["events"],
               HR = u["HR"], HR_low = u["lo"], HR_high = u["hi"], p_cox = u["p"],
               HR_adj = if (is.null(ad)) NA_real_ else ad["HR"],
               p_adj_model = if (is.null(ad)) NA_real_ else ad["p"],
               logrank_p = if (is.null(lr)) NA_real_ else 1 - pchisq(lr$chisq, 1), row.names = NULL)
  }))
  out$padj_cox <- p.adjust(out$p_cox, "BH")
  out$padj_adj_model <- p.adjust(out$p_adj_model, "BH")
  out$padj_logrank <- p.adjust(out$logrank_p, "BH")
  out$survival_direction <- ifelse(out$padj_cox < 0.05,
                                   ifelse(out$HR > 1, "High methylation = worse survival", "High methylation = better survival"),
                                   "No association")
  out %>% arrange(p_cox)
}

# ── 5. Gene level and probe level ──
tg <- cox_table(gb) %>% dplyr::rename(gene = feature)
tp <- cox_table(bt) %>% dplyr::rename(probe = feature) %>%
  dplyr::left_join(probe_gene, by = "probe") %>% dplyr::relocate(gene, .before = probe) %>% arrange(p_cox)
write.csv(tg, file.path(a, "41_meth_vs_survival_gene_clin_based.csv"),  row.names = FALSE)
write.csv(tp, file.path(a, "42_meth_vs_survival_probe_clin_based.csv"), row.names = FALSE)

cat("\n== Methylation vs survival (clinical-table cohort): GENE level ==\n")
print(as.data.frame(tg %>% dplyr::select(gene, n, events, mean_beta, HR, HR_low, HR_high, p_cox, padj_cox,
                                         HR_adj, padj_adj_model, survival_direction)), digits = 3)
cat("\n== PROBE level:", nrow(tp), "probes tested,", sum(tp$padj_cox < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")
print(as.data.frame(head(tp %>% dplyr::select(gene, probe, n, events, mean_beta, HR, p_cox, padj_cox, HR_adj, p_adj_model), 15)), digits = 3)

# ── 6. Checks against the old files and against script 61 ──
old_g <- tryCatch(read.csv(file.path(a, "9_meth_vs_survival_gene.csv")), error = function(e) NULL)
old_p <- tryCatch(read.csv(file.path(a, "10_meth_vs_survival_probe.csv")), error = function(e) NULL)
cat("\n== Old vs new (script 46 vs script 71) ==\n")
if (!is.null(old_g)) cat("Gene level: significant before =", sum(old_g$padj_cox < 0.05, na.rm = TRUE),
                         "| now =", sum(tg$padj_cox < 0.05, na.rm = TRUE), "| old n =", unique(old_g$n), "| new n =", unique(tg$n), "\n")
if (!is.null(old_p)) cat("Probe level: significant before =", sum(old_p$padj_cox < 0.05, na.rm = TRUE),
                         "| now =", sum(tp$padj_cox < 0.05, na.rm = TRUE), "\n")
cat("\ncg23932332 (DUSP10) in the new probe table:\n")
print(as.data.frame(tp %>% dplyr::filter(probe == "cg23932332") %>% dplyr::select(gene, probe, n, events, HR, p_cox, padj_cox, HR_adj, p_adj_model)), digits = 3)
f31b <- file.path(a, "31b_methylation_survival_clin_based.csv")
if (file.exists(f31b)) {
  cat("\nScript 61 (31b) for comparison; gene-level HR can differ a little because 61 uses the 17 DUSP10 probes of file 19:\n")
  print(read.csv(f31b)[c(1, 4), c("model", "HR", "p", "n", "events")], digits = 3)
}

# ── 7. Kaplan-Meier: top 3 genes + top 3 probes (median split) ──
top <- list(); for (g in head(tg$gene, 3)) top[[paste("gene", g)]] <- gb[g, ]
for (p in head(tp$probe, 3)) top[[paste("probe", p)]] <- bt[p, ]
png(file.path(fig_dir, "fig_km_methylation_top_clin_based.png"), width = 2400, height = 1600, res = 250)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (nm in names(top)) {
  x <- top[[nm]]; k <- is.finite(x); d <- clin[k, ]
  d$grp <- factor(ifelse(x[k] > median(x[k]), "High", "Low"), c("Low", "High"))
  lr <- survdiff(Surv(time, event) ~ grp, data = d)
  plot(survfit(Surv(time, event) ~ grp, data = d), col = c("#2C7BB6", "#D7301F"), lwd = 2,
       xlab = "Months", ylab = "Overall survival",
       main = sprintf("%s\nlog-rank p = %.3g", nm, 1 - pchisq(lr$chisq, 1)), cex.main = 0.9)
  legend("bottomleft", c("Low methylation", "High methylation"), col = c("#2C7BB6", "#D7301F"), lwd = 2, bty = "n", cex = 0.7)
}
dev.off()