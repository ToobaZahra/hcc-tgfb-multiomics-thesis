# 31_glmnet_signature.R
# CORRECTED VERSION:
#  1. Train/test split added -- AUC now computed on HELD-OUT samples,
#     not the same samples used to fit the model (was circular)
#  2. type.measure = "auc" set explicitly for cv.glmnet tuning
#  3. lambda.1se tried first, falls back to lambda.min if it selects
#     zero genes (same logic as 26_glmnet_survival.R)
#  4. Comparison figure: in-sample vs held-out ROC curves

library(glmnet)
library(SummarizedExperiment)
library(pROC)
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

# ------------------------------------------------------------
# Candidate genes restricted to the 7 with univariate Cox
# pval<0.05, same rationale as 26_glmnet_survival.R: keeps the
# events/samples-per-predictor ratio sane instead of throwing
# all 26 genes at a few hundred samples.
# ------------------------------------------------------------
cox_univariate <- read.csv("03_results/survival/cox_univariate.csv")
nexus_genes <- cox_univariate$gene[cox_univariate$pval < 0.05]
cat("Candidate genes (univariate pval<0.05):\n")
print(nexus_genes)

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

cat("Class balance:\n")
print(table(y))

# ------------------------------------------------------------
# Train/test split -- fit on train only, evaluate on held-out test
# ------------------------------------------------------------
set.seed(42)
n <- nrow(nexus_expr)
train_idx <- sample(seq_len(n), size = round(0.7 * n))

x_train <- nexus_expr[train_idx, ]; y_train <- y[train_idx]
x_test  <- nexus_expr[-train_idx, ]; y_test  <- y[-train_idx]

cat("Train samples:", length(y_train), " Test samples:", length(y_test), "\n")
cat("Train class balance:\n"); print(table(y_train))
cat("Test class balance:\n");  print(table(y_test))

# ------------------------------------------------------------
# glmnet logistic regression, tuned on AUC explicitly
# ------------------------------------------------------------
cv_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1,
                    nfolds = 10, type.measure = "auc")

cat("lambda.min:", cv_fit$lambda.min, "  lambda.1se:", cv_fit$lambda.1se, "\n")

lambda_used <- "lambda.1se"
coefs <- coef(cv_fit, s = "lambda.1se")
coefs_df <- data.frame(gene = rownames(coefs), coef = as.numeric(coefs))
coefs_df <- coefs_df[coefs_df$coef != 0 & coefs_df$gene != "(Intercept)", ]

if (nrow(coefs_df) == 0) {
  cat("lambda.1se selected zero genes -- falling back to lambda.min.\n",
      "NOTE: report this explicitly -- lambda.min is less conservative.\n")
  lambda_used <- "lambda.min"
  coefs <- coef(cv_fit, s = "lambda.min")
  coefs_df <- data.frame(gene = rownames(coefs), coef = as.numeric(coefs))
  coefs_df <- coefs_df[coefs_df$coef != 0 & coefs_df$gene != "(Intercept)", ]
}

coefs_df <- coefs_df[order(abs(coefs_df$coef), decreasing = TRUE), ]

cat("\nSelected genes for stage signature (fit on training set, using",
    lambda_used, "):\n")
print(coefs_df)

if (nrow(coefs_df) == 0) {
  stop("Both lambda.1se and lambda.min selected zero genes -- no signature ",
       "can be extracted from this training split. Report as a null result.")
}

write.csv(coefs_df, "03_results/signature/glmnet_stage_signature.csv", row.names = FALSE)
cat("Signature fit using:", lambda_used, "-- record this in your methods section.\n")

# ------------------------------------------------------------
# Held-out AUC -- the honest validation number
# ------------------------------------------------------------
pred_test <- predict(cv_fit, newx = x_test, s = lambda_used, type = "response")
roc_test <- roc(y_test, as.numeric(pred_test), quiet = TRUE)
auc_test <- as.numeric(auc(roc_test))
cat("\nHeld-out test AUC:", round(auc_test, 3), "on", length(y_test), "test samples\n")

# ------------------------------------------------------------
# In-sample AUC for comparison -- inflated, NOT the number to report
# ------------------------------------------------------------
pred_train <- predict(cv_fit, newx = x_train, s = lambda_used, type = "response")
roc_train <- roc(y_train, as.numeric(pred_train), quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
cat("In-sample (training) AUC:", round(auc_train, 3),
    "-- inflated, do not report as validation\n")

# ------------------------------------------------------------
# Comparison figure: in-sample vs held-out ROC curves
# ------------------------------------------------------------
png("04_figures/glmnet_stage_roc_insample_vs_heldout.png",
    width = 8, height = 6, units = "in", res = 300)
plot(roc_train, col = "#D9534F", lwd = 2,
     main = "Stage Signature: In-Sample vs Held-Out ROC")
lines(roc_test, col = "#2E7C4A", lwd = 2)
legend("bottomright",
       legend = c(paste0("In-sample (train), AUC=", round(auc_train, 3)),
                  paste0("Held-out (test), AUC=", round(auc_test, 3))),
       col = c("#D9534F", "#2E7C4A"), lwd = 2)
dev.off()

cat("\nComparison ROC figure saved to",
    "04_figures/glmnet_stage_roc_insample_vs_heldout.png\n")
cat("Report the HELD-OUT AUC (", round(auc_test, 3),
    ") as the validation result, not the in-sample number.\n")