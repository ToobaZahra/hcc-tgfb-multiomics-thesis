## 46: promoter methylation vs overall survival -- gene level (mean promoter beta) and probe level, tumors only
## Cox regression (HR per 1 SD of beta), adjusted model (age + stage), median-split log-rank, Kaplan-Meier
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(survival)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; dir.create(a, showWarnings = FALSE, recursive = TRUE); fig_dir <- "04_figures"

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Promoter probes of the 26 genes ───────────────────
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

# ── 2. Tumors + survival data (same definitions as script 45) ──
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cdm <- as.data.frame(colData(se))[tum, ]; cdm$patient <- substr(colnames(se)[tum], 1, 12)
need <- c("vital_status", "days_to_death", "days_to_last_follow_up")
if (length(setdiff(need, colnames(cdm))) > 0) stop("Missing columns: ", paste(setdiff(need, colnames(cdm)), collapse = ", "))
age_col <- intersect(c("age_at_index", "age_at_diagnosis"), colnames(cdm))[1]
age <- if (is.na(age_col)) NA_real_ else as.numeric(cdm[[age_col]])
if (!is.na(age_col) && median(age, na.rm = TRUE) > 200) age <- age / 365.25
dead <- tolower(cdm$vital_status) == "dead"
time <- as.numeric(ifelse(dead, cdm$days_to_death, cdm$days_to_last_follow_up)) / 30.44
stg <- ifelse(grepl("^Stage I{1,2}[ABC]?$", cdm$ajcc_pathologic_stage), "Early",
              ifelse(grepl("^Stage (III|IV)", cdm$ajcc_pathologic_stage), "Late", NA))
clin <- data.frame(patient = cdm$patient, time = time, event = as.integer(dead),
                   age = age, stage = factor(stg, c("Early", "Late")))
clin <- clin[!is.na(clin$time) & clin$time > 0, ]
cat("Patients with survival data:", nrow(clin), "| deaths:", sum(clin$event), "\n")

# ── 3. Beta matrices ─────────────────────────────────────
bt <- assay(se[unique(probe_gene$probe), tum]); colnames(bt) <- cdm$patient
bt <- bt[, clin$patient, drop = FALSE]
gb <- t(sapply(unique(probe_gene$gene), function(gn)
  colMeans(bt[probe_gene$probe[probe_gene$gene == gn], , drop = FALSE], na.rm = TRUE)))

# ── 4. Cox helpers ───────────────────────────────────────
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

# ── 5. Gene level and probe level ────────────────────────
tg <- cox_table(gb) %>% dplyr::rename(gene = feature)
tp <- cox_table(bt) %>% dplyr::rename(probe = feature) %>%
  dplyr::left_join(probe_gene, by = "probe") %>% dplyr::relocate(gene, .before = probe) %>% arrange(p_cox)
write.csv(tg, file.path(a, "9_meth_vs_survival_gene.csv"),  row.names = FALSE)
write.csv(tp, file.path(a, "10_meth_vs_survival_probe.csv"), row.names = FALSE)

cat("\n== Methylation vs survival: GENE level ==\n")
print(as.data.frame(tg %>% dplyr::select(gene, mean_beta, HR, HR_low, HR_high, p_cox, padj_cox,
                                         HR_adj, padj_adj_model, padj_logrank, survival_direction)), digits = 3)
cat("\n== PROBE level:", nrow(tp), "probes tested,", sum(tp$padj_cox < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")
print(as.data.frame(head(tp %>% dplyr::select(gene, probe, mean_beta, HR, p_cox, padj_cox, HR_adj, p_adj_model), 15)), digits = 3)

# ── 6. Kaplan-Meier: top 3 genes + top 3 probes ──────────
top <- list(); for (g in head(tg$gene, 3)) top[[paste("gene", g)]] <- gb[g, ]
for (p in head(tp$probe, 3)) top[[paste("probe", p)]] <- bt[p, ]
png(file.path(fig_dir, "fig_km_methylation_top.png"), width = 2400, height = 1600, res = 250)
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