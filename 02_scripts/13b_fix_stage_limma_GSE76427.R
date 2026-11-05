# 13b_fix_stage_limma_GSE76427.R
# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/13b_fix_stage_limma_GSE76427.R
# Fix: GSE76427 matrix is on linear intensity scale -> log2 it (floor at 1), re-run
# Late vs Early limma. Backs up old files (never overwrites an existing backup) and
# writes corrected files to BOTH 03_results/ and 03_results/dge/.

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
library(GEOquery)
library(limma)

geo_dir <- "data/raw/geo"

gse76427   <- readRDS(file.path(geo_dir, "GSE76427.rds"))
eset76427  <- gse76427[[1]]
pheno76427 <- pData(eset76427)

# ---- log2 transform (same as 11b) ----
m <- exprs(eset76427)
cat("Flooring", sum(m < 1, na.rm = TRUE), "values below 1 to 1\n")
m[m < 1] <- 1
mat76427 <- log2(m)

# Map BCLC to Early/Late
pheno76427$stage_group <- ifelse(pheno76427$`bclc_staging:ch1` %in% c("0", "A"), "Early",
                                 ifelse(pheno76427$`bclc_staging:ch1` %in% c("B", "C"), "Late", NA))

keep <- pheno76427$`tissue:ch1` == "primary hepatocellular carcinoma tumor" &
  !is.na(pheno76427$stage_group)

mat_filtered <- mat76427[, keep]
stage_group  <- factor(pheno76427$stage_group[keep], levels = c("Early", "Late"))
cat("Stage group distribution:\n"); print(table(stage_group))

design  <- model.matrix(~ stage_group)
fit     <- eBayes(lmFit(mat_filtered, design))
results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

feature76427 <- fData(eset76427)
results$symbol <- feature76427[rownames(results), "Gene symbol"]

# ---- sanity check ----
cat("AveExpr summary (expect roughly 4-16):\n"); print(summary(results$AveExpr))
cat("|logFC| quantiles:\n")
print(quantile(abs(results$logFC), c(.1, .5, .9, .99), na.rm = TRUE))

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

nexus_stage <- results[results$symbol %in% nexus_genes, ]
nexus_stage <- nexus_stage[order(nexus_stage$adj.P.Val), ]

# ---- back up old files, then write corrected ones to both folders ----
backup <- function(f) {
  bk <- sub("\\.csv$", "_OLD_unlogged.csv", f)
  if (file.exists(f) && !file.exists(bk)) file.copy(f, bk)
}
for (d in c("03_results", "03_results/dge")) {
  f1 <- file.path(d, "limma_GSE76427_late_vs_early.csv")
  f2 <- file.path(d, "limma_GSE76427_stage_nexus.csv")
  backup(f1); backup(f2)
  write.csv(results,     f1, row.names = TRUE)
  write.csv(nexus_stage, f2, row.names = TRUE)
}

cat("\nGSE76427 nexus genes (Late vs Early, corrected):\n")
print(nexus_stage[, c("symbol", "logFC", "adj.P.Val")])
cat("Done.\n")