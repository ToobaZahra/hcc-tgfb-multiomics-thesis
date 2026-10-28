## 55: GSE54503 replication of DUSP10 promoter methylation (paired HCC tumor vs adjacent non-tumor)
## Uses the same 17 DUSP10 promoter probes as TCGA (file 19) and checks cg23932332 specifically.
library(dplyr); library(GEOquery)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
options(timeout = 1200)
dir.create("data/raw/geo", showWarnings = FALSE, recursive = TRUE)
key <- "cg23932332"

# ---- 1. Load series matrix (saved locally so it downloads only once) ----
e <- getGEO("GSE54503", GSEMatrix = TRUE, getGPL = FALSE, destdir = "data/raw/geo")[[1]]
p <- pData(e); B <- exprs(e)
cat("Matrix:", nrow(B), "probes x", ncol(B), "samples | value range:",
    round(range(B, na.rm = TRUE), 3), "\n")
if (min(B, na.rm = TRUE) < 0 || max(B, na.rm = TRUE) > 1.5) {
  B <- 2^B / (2^B + 1); cat("Values looked like M-values: converted to beta.\n") }

# ---- 2. Tumor / non-tumor labels (CHECK the printed table against GEO) ----
chr <- grep(":ch1$", names(p), value = TRUE)
cat("\nSample annotation columns (first values):\n")
for (cc in chr) cat(" ", cc, ":", paste(head(unique(as.character(p[[cc]])), 6), collapse = " | "), "\n")
txt <- tolower(apply(p[, c("title", chr), drop = FALSE], 1, paste, collapse = " "))
p$tissue <- ifelse(grepl("non-?tumou?r|adjacent|normal|benign|surrounding", txt), "Normal",
                   ifelse(grepl("tumou?r|hcc|carcinoma", txt), "Tumor", NA))
cat("\nTissue labels:\n"); print(table(p$tissue, useNA = "ifany"))

# ---- 3. Pairing ----
# Titles are "6089 adjacent non-tumor tissue" / "6090 HCC tumor": the sample numbers differ within a pair
# and there is no patient-ID column. So each Normal is paired with the Tumor listed right after it (GSM order),
# and the pairing is then VERIFIED with the 450K rs probes (SNP genotype probes, same patient = same genotype).
p <- p[order(rownames(p)), ]
ni <- which(p$tissue == "Normal")
okpair <- ni < nrow(p) & p$tissue[pmin(ni + 1, nrow(p))] == "Tumor"
cat("\nNormal samples directly followed by a Tumor sample:", sum(okpair), "of", length(ni), "\n")
if (sum(okpair) < 20) stop("Order-based pairing failed: paste head(p[, c('title','tissue')], 20).")
nor <- rownames(p)[ni[okpair]]; tum <- rownames(p)[ni[okpair] + 1]
cat("First pairs (normal | tumor):\n"); print(data.frame(normal = p[nor[1:4], "title"], tumor = p[tum[1:4], "title"]))

rs <- grep("^rs", rownames(B), value = TRUE)
if (length(rs) >= 10) {
  same  <- sapply(seq_along(nor), function(i) cor(B[rs, nor[i]], B[rs, tum[i]], use = "complete.obs"))
  other <- sapply(seq_along(nor), function(i) cor(B[rs, nor[i]], B[rs, tum[(i %% length(tum)) + 1]], use = "complete.obs"))
  cat("Genotype-probe correlation (", length(rs), "rs probes): assumed pairs median =", round(median(same), 3),
      "| mismatched pairs median =", round(median(other), 3), "\n")
  if (median(same) < median(other) + 0.2)
    stop("Genotype check did not confirm the order-based pairing. Do not use paired tests; paste this output.")
  cat("Pairing confirmed by genotype probes.\n")
} else cat("No rs probes in the matrix: pairing could not be verified by genotype.\n")

# ---- 4. DUSP10 promoter probes (same set as TCGA) ----
probes_all <- read.csv(file.path(a, "19_DUSP10_probe_table.csv"))$probe
probes <- probes_all[probes_all %in% rownames(B)]
cat("\nDUSP10 promoter probes found:", length(probes), "of", length(probes_all),
    "| key probe", key, "present:", key %in% probes, "\n")

paired <- function(x, y) c(delta = mean(x - y, na.rm = TRUE), p = wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value)

# gene level: mean beta across promoter probes
gene_t <- colMeans(B[probes, tum, drop = FALSE], na.rm = TRUE)
gene_n <- colMeans(B[probes, nor, drop = FALSE], na.rm = TRUE)
gl <- paired(gene_t, gene_n)
gene_level <- data.frame(dataset = "GSE54503", level = "gene (promoter mean)", n_pairs = length(tum),
                         n_probes = length(probes), delta_beta = gl["delta"], p_paired_wilcox = gl["p"], row.names = NULL)

# ---- 5. Probe level, compared with TCGA and GSE56588 ----
pl <- lapply(probes, function(pr) { r <- paired(B[pr, tum], B[pr, nor])
data.frame(probe = pr, delta_GSE54503 = r["delta"], p_paired = r["p"], row.names = NULL) }) %>% bind_rows()
pl$padj_BH <- p.adjust(pl$p_paired, "BH")
f19b <- file.path(a, "19b_DUSP10_probe_tumor_vs_normal.csv"); f17b <- file.path(a, "17b_GSE56588_probe_delta.csv")
if (file.exists(f19b)) pl <- left_join(pl, read.csv(f19b) %>% select(probe, delta_TCGA = delta), by = "probe")
if (file.exists(f17b)) pl <- left_join(pl, read.csv(f17b) %>% filter(gene == "DUSP10") %>%
                                         select(probe, delta_GSE56588 = delta), by = "probe")
pl$dist_TSS <- read.csv(file.path(a, "19_DUSP10_probe_table.csv"))$dist_TSS[match(pl$probe, probes_all)]
pl <- pl %>% arrange(dist_TSS)

write.csv(gene_level, file.path(a, "25_GSE54503_DUSP10_gene_level.csv"), row.names = FALSE)
write.csv(pl, file.path(a, "25b_GSE54503_DUSP10_probe_level.csv"), row.names = FALSE)

cat("\n== DUSP10 promoter methylation, paired tumor vs adjacent (GSE54503) ==\n"); print(gene_level, digits = 3)
cat("\n== Probe level (delta = tumor minus adjacent) ==\n"); print(as.data.frame(pl), digits = 3)
if (key %in% pl$probe) { k <- pl[pl$probe == key, ]
cat("\nKey probe", key, ": GSE54503 delta =", round(k$delta_GSE54503, 3), "| paired p =", signif(k$p_paired, 3),
    "| TCGA delta =", round(k$delta_TCGA, 3), "| GSE56588 delta =", round(k$delta_GSE56588, 3), "\n") }
cat("\nNote: 'normal' here is adjacent non-tumor tissue, not healthy liver (same as TCGA normals).\n")
cat("Clinical columns (grade, stage, survival) are NOT used: check names(pData) above for any.\n")