# 11b_fix_GSE76427_log2.R
# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/11b_fix_GSE76427_log2.R
# Fix: GSE76427 matrix is on linear intensity scale -> logFC was not log2.
# Log2-transform, re-run limma, back up old outputs, overwrite with corrected ones.

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
library(GEOquery)
library(limma)

geo_dir <- "data/raw/geo"

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

gse76427  <- readRDS(file.path(geo_dir, "GSE76427.rds"))
eset76427 <- gse76427[[1]]
m <- exprs(eset76427)

cat("Raw matrix summary:\n"); print(summary(as.vector(m)))
n_nonpos <- sum(m <= 0, na.rm = TRUE)
cat("Values <= 0:", n_nonpos, "\n")

# log2 needs positive values. Only ~0.07% of values are <= 0 (background-
# subtracted noise), so floor them at 1 instead of shifting the whole matrix.
cat("Flooring", sum(m < 1, na.rm = TRUE), "values below 1 to 1\n")
m[m < 1] <- 1
mat <- log2(m)
cat("Log2 matrix summary (expect roughly 4-16):\n"); print(summary(as.vector(mat)))

# ---- limma (same design as script 11) ----
pheno <- pData(eset76427)
tissue <- ifelse(pheno[["tissue:ch1"]] == "primary hepatocellular carcinoma tumor", "Tumor",
                 ifelse(pheno[["tissue:ch1"]] == "adjacent non-tumor liver tissue", "Normal", NA))
keep <- !is.na(tissue)
tissue_group <- factor(tissue[keep], levels = c("Normal", "Tumor"))
print(table(tissue_group))

design <- model.matrix(~ tissue_group)
fit <- eBayes(lmFit(mat[, keep], design))
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

fd <- fData(eset76427)
symcol <- intersect(c("Gene symbol", "Gene Symbol"), colnames(fd))[1]
res$symbol <- fd[rownames(res), symcol]

# ---- sanity check ----
cat("AveExpr summary (expect roughly 4-16):\n"); print(summary(res$AveExpr))
cat("|logFC| quantiles (expect median well under 1):\n")
print(quantile(abs(res$logFC), c(.1, .5, .9, .99), na.rm = TRUE))

# ---- back up old outputs, then overwrite ----
out_full  <- "03_results/limma_GSE76427_tumor_vs_normal.csv"
out_nexus <- "03_results/limma_GSE76427_nexus_genes.csv"
for (f in c(out_full, out_nexus))
{ bk <- sub("\\.csv$", "_OLD_unlogged.csv", f)
if (file.exists(f) && !file.exists(bk)) file.copy(f, bk) }  # never overwrite an existing backup

nexus <- res[res$symbol %in% nexus_genes, ]
nexus <- nexus[order(nexus$adj.P.Val), ]
write.csv(res,   out_full,  row.names = TRUE)
write.csv(nexus, out_nexus, row.names = TRUE)

cat("\nGSE76427 nexus genes (corrected):\n")
print(nexus[, c("symbol", "logFC", "adj.P.Val")])
cat("Done.\n")