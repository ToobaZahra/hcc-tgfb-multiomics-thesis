# 11_limma_GSE36376_GSE76427.R
library(GEOquery)
library(limma)

geo_dir <- "data/raw/geo"

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── Helper function ───────────────────────────────────────
run_limma <- function(mat, pheno, feature_data, tissue_col,
                      tumor_label, normal_label, gse_name) {
  pheno$tissue_group <- ifelse(pheno[[tissue_col]] == tumor_label, "Tumor",
                               ifelse(pheno[[tissue_col]] == normal_label, "Normal", NA))
  keep <- !is.na(pheno$tissue_group)
  mat_filtered <- mat[, keep]
  tissue_group <- factor(pheno$tissue_group[keep], levels = c("Normal", "Tumor"))
  
  cat(gse_name, "sample counts:\n")
  print(table(tissue_group))
  
  design <- model.matrix(~ tissue_group)
  fit <- lmFit(mat_filtered, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  
  # Map gene symbols
  if("Gene symbol" %in% colnames(feature_data)) {
    results$symbol <- feature_data[rownames(results), "Gene symbol"]
  } else if("Gene Symbol" %in% colnames(feature_data)) {
    results$symbol <- feature_data[rownames(results), "Gene Symbol"]
  }
  
  return(results)
}

# ── GSE36376 ──────────────────────────────────────────────
gse36376 <- readRDS(file.path(geo_dir, "GSE36376.rds"))
mat36376 <- readRDS(file.path(geo_dir, "GSE36376_norm.rds"))
eset36376 <- gse36376[[1]]

res36376 <- run_limma(
  mat = mat36376,
  pheno = pData(eset36376),
  feature_data = fData(eset36376),
  tissue_col = "tissue:ch1",
  tumor_label = "liver tumor",
  normal_label = "adjacent non-tumor liver",
  gse_name = "GSE36376"
)

nexus36376 <- res36376[res36376$symbol %in% nexus_genes, ]
nexus36376 <- nexus36376[order(nexus36376$adj.P.Val), ]

write.csv(res36376, "03_results/limma_GSE36376_tumor_vs_normal.csv", row.names = TRUE)
write.csv(nexus36376, "03_results/limma_GSE36376_nexus_genes.csv", row.names = TRUE)

cat("\nGSE36376 nexus genes:\n")
print(nexus36376[, c("symbol", "logFC", "adj.P.Val")])

# ── GSE76427 ──────────────────────────────────────────────
# Note: GSE76427 is already RSN-normalized by authors (lumi package)
# Using raw expression matrix directly — no additional normalization needed
gse76427 <- readRDS(file.path(geo_dir, "GSE76427.rds"))
mat76427 <- exprs(gse76427[[1]])  # use raw RSN-normalized data directly
eset76427 <- gse76427[[1]]
res76427 <- run_limma(
  mat = mat76427,
  pheno = pData(eset76427),
  feature_data = fData(eset76427),
  tissue_col = "tissue:ch1",
  tumor_label = "primary hepatocellular carcinoma tumor",
  normal_label = "adjacent non-tumor liver tissue",
  gse_name = "GSE76427"
)
nexus76427 <- res76427[res76427$symbol %in% nexus_genes, ]
nexus76427 <- nexus76427[order(nexus76427$adj.P.Val), ]
write.csv(res76427, "03_results/limma_GSE76427_tumor_vs_normal.csv", row.names = TRUE)
write.csv(nexus76427, "03_results/limma_GSE76427_nexus_genes.csv", row.names = TRUE)
cat("\nGSE76427 nexus genes:\n")
print(nexus76427[, c("symbol", "logFC", "adj.P.Val")])
cat("Done.\n")