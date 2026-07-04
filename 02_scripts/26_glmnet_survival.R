# 26_glmnet_survival.R
library(glmnet)
library(survival)
library(SummarizedExperiment)
library(dplyr)

# Load data
clin_tumor <- readRDS("03_results/survival_clinical_matched.rds")
vsd_tumor <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds")
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# Build expression matrix for nexus genes only
expr_mat <- assay(vsd_tumor)
gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

# Extract nexus gene expression
nexus_expr <- matrix(NA, nrow = ncol(vsd_tumor), ncol = nrow(gene_map))
colnames(nexus_expr) <- gene_map$symbol
rownames(nexus_expr) <- colnames(vsd_tumor)

for (i in 1:nrow(gene_map)) {
  ens <- gene_map$ensembl_id[i]
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) > 0) {
    nexus_expr[, i] <- expr_mat[row_match[1], ]
  }
}

# Remove columns with all NA
nexus_expr <- nexus_expr[, colSums(is.na(nexus_expr)) == 0]
cat("Expression matrix:", nrow(nexus_expr), "samples x", ncol(nexus_expr), "genes\n")

# Survival object
surv_obj <- Surv(clin_tumor$os_time, clin_tumor$os_event)

# glmnet cox
set.seed(42)
cv_fit <- cv.glmnet(nexus_expr, surv_obj, 
                    family = "cox", 
                    alpha = 1,
                    nfolds = 10)

cat("Best lambda:", cv_fit$lambda.min, "\n")

# Extract coefficients
coefs <- coef(cv_fit, s = "lambda.min")
coefs_df <- data.frame(
  gene = rownames(coefs),
  coefficient = as.numeric(coefs)
)
coefs_df <- coefs_df[coefs_df$coefficient != 0, ]
coefs_df <- coefs_df[order(abs(coefs_df$coefficient), decreasing = TRUE), ]

cat("\nSelected genes in survival signature:\n")
print(coefs_df)

write.csv(coefs_df, "03_results/glmnet_survival_signature.csv", row.names = FALSE)