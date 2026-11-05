## 67: Kaplan-Meier survival plots (median split) for DUSP10, SMAD2, MYC, MAPK8, HDAC11
##   Expression OS : TCGA-LIHC (VST), GSE14520 (log2), GSE76427 (log2)
##   Expression RFS: GSE14520, GSE76427
##   Methylation OS: TCGA-LIHC only (gene-level mean beta, and strongest probe)
## Each panel: KM curves, log-rank p, Cox HR (High vs Low) with 95% CI, n and events.
## GEO genes with several probes use the mean across the gene's probes.
## Saves PNGs to 04_figures and the numbers to 03_results/association/37_survival_5genes_stats.csv
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(GEOquery); library(survival); library(ggplot2); library(dplyr) })

# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/67_survival_km_5genes.R
setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
fig_dir <- "04_figures"; a <- "03_results/association"; g <- "data/raw/geo/"
dir.create(fig_dir, showWarnings = FALSE)
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")
cols  <- c(Low = "#2C7BB6", High = "#D7301F")

# ---------- helpers ----------
pickc <- function(df, cands, what) {
  h <- intersect(cands, names(df))
  if (length(h) == 0) stop("No ", what, " column. Columns found: ", paste(names(df), collapse = ", "))
  h[1] }
find_file <- function(pat) {
  f <- list.files(c("data/raw", "03_results"), pattern = pat, recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) stop("File not found: ", pat); f[1] }

km_one <- function(d, row, col) {
  d <- d[is.finite(d$time) & is.finite(d$event) & is.finite(d$x) & d$time > 0, ]
  if (nrow(d) < 20 || sum(d$event) < 5) { cat("Skipped (too few samples/events):", row, col, "\n"); return(NULL) }
  d$group <- factor(ifelse(d$x > median(d$x), "High", "Low"), c("Low", "High"))
  if (length(unique(d$group)) < 2) return(NULL)
  sf <- survfit(Surv(time, event) ~ group, data = d)
  km <- data.frame(time = sf$time, surv = sf$surv, nc = sf$n.censor,
                   group = rep(sub("group=", "", names(sf$strata)), sf$strata))
  km <- rbind(data.frame(time = 0, surv = 1, nc = 0, group = c("Low", "High")), km)
  km$row <- row; km$col <- col
  lr <- survdiff(Surv(time, event) ~ group, data = d); p <- 1 - pchisq(lr$chisq, 1)
  ci <- summary(coxph(Surv(time, event) ~ group, data = d))$conf.int[1, ]
  list(km = km, st = data.frame(row = row, col = col, n = nrow(d), events = sum(d$event),
                                hr = ci[1], lo = ci[3], hi = ci[4], p_logrank = p)) }

plot_km <- function(res, row_levels, col_levels, title, file, height, width = 10) {
  res <- Filter(Negate(is.null), res); if (length(res) == 0) { cat("Nothing to plot for", file, "\n"); return(NULL) }
  km <- bind_rows(lapply(res, `[[`, "km")); st <- bind_rows(lapply(res, `[[`, "st"))
  km$row <- factor(km$row, row_levels); km$col <- factor(km$col, col_levels); km$group <- factor(km$group, c("Low", "High"))
  st$row <- factor(st$row, row_levels); st$col <- factor(st$col, col_levels)
  st$txt <- sprintf("HR %.2f (%.2f-%.2f)\nlog-rank p = %.2g\nn = %d, events = %d", st$hr, st$lo, st$hi, st$p_logrank, st$n, st$events)
  p <- ggplot(km, aes(time, surv, color = group)) +
    geom_step(linewidth = 0.7) +
    geom_point(data = km[km$nc > 0, ], shape = 3, size = 1) +
    geom_text(data = st, aes(x = -Inf, y = -Inf, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = -0.15, size = 2.6) +
    facet_grid(row ~ col, scales = "free_x") +
    scale_color_manual(values = cols, name = "Median split") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = title, subtitle = "Kaplan-Meier; log-rank (raw p); HR = High vs Low (Cox)", x = "Time (months)", y = "Survival probability") +
    theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(fig_dir, file), p, width = width, height = height, dpi = 300, bg = "white")
  st }

# ---------- TCGA clinical (OS) ----------
clin <- as.data.frame(readRDS(find_file("^tcga_lihc_clinical\\.rds$")))
idc  <- pickc(clin, c("submitter_id", "bcr_patient_barcode", "patient", "case_submitter_id"), "patient ID")
vc   <- pickc(clin, c("vital_status", "Vital_status"), "vital status")
dtd  <- pickc(clin, c("days_to_death", "Days_to_death"), "days to death")
dlf  <- pickc(clin, c("days_to_last_follow_up", "days_to_last_followup", "Days_to_last_follow_up"), "days to last follow-up")
dead <- tolower(as.character(clin[[vc]])) %in% c("dead", "deceased", "1")
tcga_os <- data.frame(patient = substr(as.character(clin[[idc]]), 1, 12),
                      time = as.numeric(ifelse(dead, clin[[dtd]], clin[[dlf]])) / 30.4375,
                      event = as.numeric(dead))
tcga_os <- tcga_os[!duplicated(tcga_os$patient), ]
cat("TCGA clinical:", nrow(tcga_os), "patients,", sum(tcga_os$event, na.rm = TRUE), "deaths\n")

# ---------- TCGA expression ----------
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
pat_e <- substr(colnames(expr)[et], 1, 12)
ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw <- match(ens, sub("\\..*$", "", rownames(expr)))
if (any(is.na(rw))) stop("Not found in the TCGA expression matrix: ", paste(genes[is.na(rw)], collapse = ", "))

# ---------- GEO expression + survival ----------
prep_annot <- function(f) { tab <- Table(getGEO(filename = f))
sc <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
list(ids = as.character(tab$ID), sy = lapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), trimws)) }
gene_mean <- function(ann, mat, gene) {
  idx <- which(vapply(ann$sy, function(x) gene %in% x, logical(1)))
  hit <- intersect(ann$ids[idx], rownames(mat)); if (length(hit) == 0) return(NULL)
  m <- mat[hit, , drop = FALSE]
  if (max(mat, na.rm = TRUE) > 100) { m[m < 1] <- 1; m <- log2(m) }   # linear scale -> floor at 1, log2
  colMeans(m, na.rm = TRUE) }

m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
tis <- tolower(trimws(cl$`Tissue Type`))
s1 <- data.frame(id = as.character(cl$Affy_GSM), tumor = !grepl("non", tis) & grepl("tumor", tis),
                 os_t = as.numeric(cl$`Survival months`), os_e = as.numeric(cl$`Survival status`),
                 rfs_t = as.numeric(cl$`Recurr months`), rfs_e = as.numeric(cl$`Recurr status`))
s1 <- s1[s1$tumor & s1$id %in% colnames(m1) & !duplicated(s1$id), ]
cat("GSE14520 tumors with data:", nrow(s1), "| OS event values:", paste(names(table(s1$os_e)), collapse = "/"),
    "| RFS event values:", paste(names(table(s1$rfs_e)), collapse = "/"), "\n")
ann1 <- prep_annot(paste0(g, "GPL3921.annot.gz"))

m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]]); p2 <- p2[match(colnames(m2), rownames(p2)), ]
tis2 <- tolower(as.character(p2$`tissue:ch1`))
s2 <- data.frame(id = colnames(m2), tumor = !grepl("non", tis2) & grepl("tumor", tis2),
                 os_t = as.numeric(p2$`duryears_os:ch1`) * 12, os_e = as.numeric(p2$`event_os:ch1`),
                 rfs_t = as.numeric(p2$`duryears_rfs:ch1`) * 12, rfs_e = as.numeric(p2$`event_rfs:ch1`))
s2 <- s2[s2$tumor & !duplicated(s2$id), ]
cat("GSE76427 tumors with data:", nrow(s2), "| OS event values:", paste(names(table(s2$os_e)), collapse = "/"),
    "| RFS event values:", paste(names(table(s2$rfs_e)), collapse = "/"), "\n")
ann2 <- prep_annot(paste0(g, "GPL10558.annot.gz"))

geo_d <- function(s, m, ann, G, tcol, ecol) {
  v <- gene_mean(ann, m, G); if (is.null(v)) { cat(G, "not found on the array\n"); return(NULL) }
  data.frame(time = s[[tcol]], event = s[[ecol]], x = unname(v[s$id])) }

# ---------- Expression figures ----------
os_res <- unlist(lapply(genes, function(G) {
  i <- match(G, genes)
  dT <- merge(data.frame(patient = pat_e, x = as.numeric(expr[rw[i], et])), tcga_os, by = "patient")
  list(km_one(dT, G, "TCGA-LIHC"),
       { d <- geo_d(s1, m1, ann1, G, "os_t", "os_e"); if (is.null(d)) NULL else km_one(d, G, "GSE14520") },
       { d <- geo_d(s2, m2, ann2, G, "os_t", "os_e"); if (is.null(d)) NULL else km_one(d, G, "GSE76427") })
}), recursive = FALSE)
st_os <- plot_km(os_res, genes, c("TCGA-LIHC", "GSE14520", "GSE76427"),
                 "Overall survival by expression (median split)", "fig_5genes_survival_expression_OS.png", 14)

rfs_res <- unlist(lapply(genes, function(G) list(
  { d <- geo_d(s1, m1, ann1, G, "rfs_t", "rfs_e"); if (is.null(d)) NULL else km_one(d, G, "GSE14520") },
  { d <- geo_d(s2, m2, ann2, G, "rfs_t", "rfs_e"); if (is.null(d)) NULL else km_one(d, G, "GSE76427") })), recursive = FALSE)
st_rfs <- plot_km(rfs_res, genes, c("GSE14520", "GSE76427"),
                  "Recurrence-free survival by expression (median split)", "fig_5genes_survival_expression_RFS.png", 14, 7)

# ---------- Methylation (TCGA OS) ----------
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds"); beta <- assay(se)
tum <- which(substr(colnames(se), 14, 15) == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
pat_m <- substr(colnames(se)[tum], 1, 12)
pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))
gprobes <- list(); keyp <- character(0)
for (G in genes) {
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No methylation probes for", G, ": skipped.\n"); next }
  gprobes[[G]] <- pg$probe
  keyp[G] <- if (G == "DUSP10" && "cg23932332" %in% pg$probe) "cg23932332" else pg$probe[which.max(abs(pg$rho))]
}
mg <- names(gprobes)
meth_d <- function(probes) {
  b <- if (length(probes) == 1) as.numeric(beta[probes, tum]) else colMeans(beta[probes, tum, drop = FALSE], na.rm = TRUE)
  merge(data.frame(patient = pat_m, x = b), tcga_os, by = "patient") }

mg_res <- lapply(mg, function(G) km_one(meth_d(gprobes[[G]]), G, "TCGA-LIHC"))
st_mg <- plot_km(mg_res, mg, "TCGA-LIHC", "Overall survival by methylation, gene-level mean beta (median split)",
                 "fig_5genes_survival_methylation_gene_level_TCGA.png", 14, 5)
mp_res <- lapply(mg, function(G) km_one(meth_d(keyp[[G]]), paste0(G, "\n", keyp[[G]]), "TCGA-LIHC"))
st_mp <- plot_km(mp_res, paste0(mg, "\n", keyp[mg]), "TCGA-LIHC", "Overall survival by methylation, strongest probe (median split)",
                 "fig_5genes_survival_methylation_key_probe_TCGA.png", 14, 5)

# ---------- Table ----------
tab <- bind_rows(
  if (!is.null(st_os))  st_os  %>% mutate(measure = "Expression OS"),
  if (!is.null(st_rfs)) st_rfs %>% mutate(measure = "Expression RFS"),
  if (!is.null(st_mg))  st_mg  %>% mutate(measure = "Methylation OS (gene-level)"),
  if (!is.null(st_mp))  st_mp  %>% mutate(measure = "Methylation OS (strongest probe)", row = sub("\n.*", "", as.character(row))))
write.csv(tab, file.path(a, "37_survival_5genes_stats.csv"), row.names = FALSE)
cat("\n== Survival results (median split) ==\n"); print(as.data.frame(tab), digits = 3)
cat("\nSaved 4 figures to", fig_dir, "and 37_survival_5genes_stats.csv\n")
cat("Raw log-rank p-values, univariate, no adjustment for stage or age. Methylation survival is TCGA only.\n")