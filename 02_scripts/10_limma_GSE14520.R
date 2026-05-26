# 10_limma_GSE14520.R
library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)

geo_dir <- "data/raw/geo"

# Load preprocessed GEO data
gse14520 <- readRDS(file.path(geo_dir, "GSE14520.rds"))
mat_norm <- readRDS(file.path(geo_dir, "GSE14520_norm.rds"))

# Get phenotype data
eset <- gse14520[[1]]
pheno <- pData(eset)
feature_data <- fData(eset)

# Use Tissue:ch1 for tumor/normal label
pheno$tissue_group <- ifelse(pheno$`Tissue:ch1` == "Liver Tumor Tissue", "Tumor",
                             ifelse(pheno$`Tissue:ch1` == "Liver Non-Tumor Tissue", "Normal", NA))

# Remove samples with no label
keep <- !is.na(pheno$tissue_group)
mat_filtered <- mat_norm[, keep]
tissue_group <- factor(pheno$tissue_group[keep], levels = c("Normal", "Tumor"))

cat("Sample counts:\n")
print(table(tissue_group))

# Design matrix
design <- model.matrix(~ tissue_group)

# Fit limma model
fit <- lmFit(mat_filtered, design)
fit <- eBayes(fit)

# Get results
results_limma <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Map probe IDs to gene symbols
results_limma$symbol <- feature_data[rownames(results_limma), "Gene symbol"]

# Filter nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

nexus_results <- results_limma[results_limma$symbol %in% nexus_genes, ]
nexus_results <- nexus_results[order(nexus_results$adj.P.Val), ]

# Save
write.csv(results_limma, "03_results/limma_GSE14520_tumor_vs_normal.csv", row.names = TRUE)
write.csv(nexus_results, "03_results/limma_GSE14520_nexus_genes.csv", row.names = TRUE)

cat("\nNexus genes in GSE14520:\n")
print(nexus_results[, c("symbol", "logFC", "adj.P.Val")])
cat("Done.\n")