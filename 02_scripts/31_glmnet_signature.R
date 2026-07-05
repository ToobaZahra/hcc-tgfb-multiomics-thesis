# 31_glmnet_signature.R
library(glmnet)
library(SummarizedExperiment)
library(dplyr)

vsd <- readRDS("data/raw/tcga_lihc_vsd.rds")
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")
res_df <- readRDS("03_results/dge/dge_tcga_tumor_vs_normal.rds")

# Add stage group
clin$stage_group <- ifelse(
  clin$ajcc_pathologic_stage %in% c("Stage I", "Stage II"), "Early",
  ifelse(clin$ajcc_pathologic_stage %in%
           c("Stage III","Stage IIIA","Stage IIIB","Stage IV"), "Late", NA))

# Match samples
vsd$patient <- substr(colnames(vsd), 1, 12)
tumor_keep <- vsd$sample_type == "Primary Tumor"
vsd_tumor <- vsd[, tumor_keep]
common <- intersect(clin$submitter_id, vsd_tumor$patient)
clin_matched <- clin[clin$submitter_id %in% common, ]
vsd_matched <- vsd_tumor[, vsd_tumor$patient %in% common]
clin_matched <- clin_matched[match(vsd_matched$patient, clin_matched$submitter_id), ]

# Keep only staged samples
keep <- !is.na(clin_matched$stage_group)
clin_matched <- clin_matched[keep, ]
vsd_matched <- vsd_matched[, keep]
cat("Staged samples:", ncol(vsd_matched), "\n")
cat("Stage distribution:\n")
print(table(clin_matched$stage_group))








# Build expression matrix for nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

expr_mat <- assay(vsd_matched)
nexus_expr <- matrix(NA, nrow = ncol(vsd_matched), ncol = nrow(gene_map))
colnames(nexus_expr) <- gene_map$symbol

for (i in 1:nrow(gene_map)) {
  ens <- gene_map$ensembl_id[i]
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) > 0) nexus_expr[, i] <- expr_mat[row_match[1], ]
}

nexus_expr <- nexus_expr[, colSums(is.na(nexus_expr)) == 0]
y <- factor(clin_matched$stage_group, levels = c("Early", "Late"))

# glmnet logistic regression
set.seed(42)
cv_fit <- cv.glmnet(nexus_expr, y, family = "binomial", alpha = 1, nfolds = 10)

coefs <- coef(cv_fit, s = "lambda.min")
coefs_df <- data.frame(gene = rownames(coefs), coef = as.numeric(coefs))
coefs_df <- coefs_df[coefs_df$coef != 0 & coefs_df$gene != "(Intercept)", ]
coefs_df <- coefs_df[order(abs(coefs_df$coef), decreasing = TRUE), ]

cat("Selected genes for stage signature:\n")
print(coefs_df)

saveRDS(cv_fit, "03_results/signature/glmnet_fit.rds")
write.csv(coefs_df, "03_results/signature/glmnet_stage_signature.csv", row.names = FALSE)