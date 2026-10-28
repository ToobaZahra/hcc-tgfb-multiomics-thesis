## 57: DUSP10 methylation split by genomic region
## Question: is the tumor-vs-normal signal in the main-TSS CpG cluster or only in the distal probes?
## Sets: main-TSS cluster (all / unmasked), distal probes, and cg23932332 alone.
## Cohorts: TCGA-LIHC (tumor vs normal + correlation with expression) and GSE54503 (paired), if cached locally.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(org.Hs.eg.db) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
G <- "DUSP10"; key <- "cg23932332"
tss <- 221742090   # hg38 chr1; midpoint of cg27632268 (+5) and cg08708209 (-8). CHECK strand and TSS in a genome browser.
win <- 1500        # probes within +/- 1.5 kb of that TSS = "main-TSS cluster"

# ---- 1. Probe annotation from the TCGA object (hg38 coordinates + quality masks) ----
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
rd <- as.data.frame(rowData(se))
pt <- read.csv(file.path(a, "19_DUSP10_probe_table.csv"))
x <- rd[pt$probe, ]
ann <- data.frame(probe = pt$probe, dist_TSS_file = pt$dist_TSS, mean_beta_tumor = pt$mean_beta,
                  pos_hg38 = x$beg_A, probeCpGcnt = x$probeCpGcnt, context35 = x$context35,
                  masked = (x$MASK_general %in% TRUE) | (x$MASK_snp5_common %in% TRUE))
ann$dist_main_TSS <- ann$pos_hg38 - tss
ann$region <- ifelse(abs(ann$dist_main_TSS) <= win, "main-TSS cluster", "distal")

sets <- list(
  `main-TSS cluster (all)`      = ann$probe[ann$region == "main-TSS cluster"],
  `main-TSS cluster (unmasked)` = ann$probe[ann$region == "main-TSS cluster" & !ann$masked],
  `distal probes`               = ann$probe[ann$region == "distal"],
  `cg23932332 alone`            = key)
cat("Probes per set:\n"); print(sapply(sets, length))
cat("Distal probes:", paste(sets[["distal probes"]], collapse = ", "), "\n")
cat("Masked probes:", paste(ann$probe[ann$masked], collapse = ", "), "\n")

# ---- 2. TCGA-LIHC: tumor vs normal and methylation vs expression, per set ----
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
M <- assay(se[ann$probe, ])
pat <- substr(colnames(se)[tum], 1, 12)

vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ex <- expr[row, et]; names(ex) <- substr(colnames(expr)[et], 1, 12)

tcga <- lapply(names(sets), function(s) {
  gl <- colMeans(M[sets[[s]], , drop = FALSE], na.rm = TRUE)
  vt <- unname(gl[tum]); vn <- unname(gl[nor]); y <- unname(ex[pat])
  k <- is.finite(vt) & is.finite(y)
  ct <- suppressWarnings(cor.test(vt[k], y[k], method = "spearman", exact = FALSE))
  data.frame(cohort = "TCGA-LIHC", set = s, n_probes = length(sets[[s]]),
             delta_beta = mean(vt, na.rm = TRUE) - mean(vn, na.rm = TRUE),
             p_tumor_vs_normal = wilcox.test(vt, vn)$p.value,
             rho_expr = unname(ct$estimate), p_rho = ct$p.value, row.names = NULL) }) %>% bind_rows()

# ---- 3. GSE54503 (paired), if the series matrix is cached from script 55 ----
res <- tcga
f54 <- "data/raw/geo/GSE54503_series_matrix.txt.gz"
if (file.exists(f54)) {
  suppressPackageStartupMessages(library(GEOquery))
  e <- getGEO(filename = f54, getGPL = FALSE); p <- pData(e); B <- exprs(e)
  chr <- grep(":ch1$", names(p), value = TRUE)
  txt <- tolower(apply(p[, c("title", chr), drop = FALSE], 1, paste, collapse = " "))
  p$tissue <- ifelse(grepl("non-?tumou?r|adjacent|normal|benign|surrounding", txt), "Normal",
                     ifelse(grepl("tumou?r|hcc|carcinoma", txt), "Tumor", NA))
  p <- p[order(rownames(p)), ]
  ni <- which(p$tissue == "Normal")
  ok <- ni < nrow(p) & p$tissue[pmin(ni + 1, nrow(p))] == "Tumor"
  nor_g <- rownames(p)[ni[ok]]; tum_g <- rownames(p)[ni[ok] + 1]
  rs <- grep("^rs", rownames(B), value = TRUE)
  same <- sapply(seq_along(nor_g), function(i) cor(B[rs, nor_g[i]], B[rs, tum_g[i]], use = "complete.obs"))
  other <- sapply(seq_along(nor_g), function(i) cor(B[rs, nor_g[i]], B[rs, tum_g[(i %% length(tum_g)) + 1]], use = "complete.obs"))
  if (median(same) > median(other) + 0.2) {
    cat("\nGSE54503: pairing confirmed by genotype probes;", length(nor_g), "pairs.\n")
    g54 <- lapply(names(sets), function(s) {
      pr <- intersect(sets[[s]], rownames(B))
      gt <- colMeans(B[pr, tum_g, drop = FALSE], na.rm = TRUE); gn <- colMeans(B[pr, nor_g, drop = FALSE], na.rm = TRUE)
      data.frame(cohort = "GSE54503 (paired)", set = s, n_probes = length(pr), delta_beta = mean(gt - gn, na.rm = TRUE),
                 p_tumor_vs_normal = wilcox.test(gt, gn, paired = TRUE, exact = FALSE)$p.value,
                 rho_expr = NA, p_rho = NA, row.names = NULL) }) %>% bind_rows()
    res <- bind_rows(res, g54)
  } else cat("\nGSE54503 skipped: genotype check did not confirm the pairing.\n")
} else cat("\nGSE54503 series matrix not found in data/raw/geo: GSE54503 rows skipped.\n")

write.csv(res, file.path(a, "27_DUSP10_methylation_by_region.csv"), row.names = FALSE)
cat("\n== DUSP10 methylation by genomic region ==\n"); print(res, digits = 3)

# ---- 4. Probe annotation with coordinates, quality flags and deltas in all cohorts ----
f23b <- file.path(a, "23b_DUSP10_methylation_probe_level.csv"); f25b <- file.path(a, "25b_GSE54503_DUSP10_probe_level.csv")
out <- ann
if (file.exists(f23b)) out <- left_join(out, read.csv(f23b) %>% dplyr::select(probe, delta_TCGA, delta_GSE56588), by = "probe")
if (file.exists(f25b)) out <- left_join(out, read.csv(f25b) %>% dplyr::select(probe, delta_GSE54503, p_paired_GSE54503 = p_paired), by = "probe")
out <- out %>% arrange(pos_hg38)
write.csv(out, file.path(a, "26c_DUSP10_probe_annotation_coords.csv"), row.names = FALSE)
cat("\n== Probes ordered by genomic position ==\n"); print(out, digits = 3)
cat("\nNotes: TSS position (", tss, ") is inferred from the probes nearest dist_TSS = 0; confirm in a genome browser.\n", sep = "")
cat("Raw p-values. 'distal' = more than ", win, " bp from that TSS. GSE56588 is not included (needs its large beta file).\n", sep = "")