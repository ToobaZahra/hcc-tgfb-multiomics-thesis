# 13_stage_limma_GSE76427.R
library(GEOquery)
library(limma)

geo_dir <- "data/raw/geo"

gse76427 <- readRDS(file.path(geo_dir, "GSE76427.rds"))
mat76427 <- exprs(gse76427[[1]])
eset76427 <- gse76427[[1]]
pheno76427 <- pData(eset76427)

# Check BCLC staging
cat("BCLC staging values:\n")
print(table(pheno76427$`bclc_staging:ch1`))
cat("\nTissue values:\n")
print(table(pheno76427$`tissue:ch1`))

# Map BCLC to Early/Late
pheno76427$stage_group <- ifelse(pheno76427$`bclc_staging:ch1` %in% c("0", "A"), "Early",
                                 ifelse(pheno76427$`bclc_staging:ch1` %in% c("B", "C"), "Late", NA))

# Tumor samples only with stage info
keep <- pheno76427$`tissue:ch1` == "primary hepatocellular carcinoma tumor" & 
  !is.na(pheno76427$stage_group)

mat_filtered <- mat76427[, keep]
stage_group <- factor(pheno76427$stage_group[keep], levels = c("Early", "Late"))

cat("Stage group distribution:\n")
print(table(stage_group))

# Design matrix
design <- model.matrix(~ stage_group)
fit <- lmFit(mat_filtered, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Map gene symbols
feature76427 <- fData(eset76427)
results$symbol <- feature76427[rownames(results), "Gene symbol"]

# Filter nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

nexus_stage <- results[results$symbol %in% nexus_genes, ]
nexus_stage <- nexus_stage[order(nexus_stage$adj.P.Val), ]

write.csv(results, "03_results/limma_GSE76427_late_vs_early.csv", row.names = TRUE)
write.csv(nexus_stage, "03_results/limma_GSE76427_stage_nexus.csv", row.names = TRUE)

cat("\nGSE76427 nexus genes (Late vs Early):\n")
print(nexus_stage[, c("symbol", "logFC", "adj.P.Val")])
cat("Done.\n")