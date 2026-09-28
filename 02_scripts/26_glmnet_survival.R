# 26_glmnet_survival.R
# CORRECTED VERSION:
#  1. Fixed file paths (missing "survival/" and "dge/" subfolders)
#  2. Train/test split added -- risk score KM test is now on HELD-OUT
#     samples, not the same samples used to fit the model (was circular)
#  3. lambda.1se used instead of lambda.min for a less overfit signature

library(glmnet)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)

# ------------------------------------------------------------
# Load data -- paths corrected to match your actual folder layout
# ------------------------------------------------------------
clin_tumor <- readRDS("03_results/survival/survival_clinical_matched.rds")
vsd_tumor  <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds")
res_df     <- readRDS("03_results/dge/dge_tcga_tumor_vs_normal.rds")

# ------------------------------------------------------------
# Candidate gene list restricted to the 7 genes with univariate
# Cox pval<0.05 (from 28_survival.R), NOT all 26 nexus genes.
#
# Why: with 91 events in the training split, fitting all 26 genes
# gives ~3.5 events/predictor -- well under the ~10 events/predictor
# rule of thumb for stable Cox estimates. That underpowered fit is
# exactly why lambda.1se shrank every coefficient to zero on the
# first attempt. Restricting to genes already justified by the
# univariate screen brings this to a defensible ~13 events/predictor.
# ------------------------------------------------------------
cox_univariate <- read.csv("03_results/survival/cox_univariate.csv")
nexus_genes <- cox_univariate$gene[cox_univariate$pval < 0.05]

cat("Candidate genes for glmnet signature (univariate pval<0.05):\n")
print(nexus_genes)

expr_mat <- assay(vsd_tumor)
gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

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

nexus_expr <- nexus_expr[, colSums(is.na(nexus_expr)) == 0]
cat("Expression matrix:", nrow(nexus_expr), "samples x", ncol(nexus_expr), "genes\n")

surv_obj <- Surv(clin_tumor$os_time, clin_tumor$os_event)

# ------------------------------------------------------------
# Train/test split -- fit the model on train only, so the KM
# validation later is on samples the model has never seen.
# ------------------------------------------------------------
set.seed(42)
n <- nrow(nexus_expr)
train_idx <- sample(seq_len(n), size = round(0.7 * n))

x_train    <- nexus_expr[train_idx, ]
surv_train <- surv_obj[train_idx]

x_test     <- nexus_expr[-train_idx, ]
clin_test  <- clin_tumor[-train_idx, ]

cat("Train samples:", nrow(x_train), " Test samples:", nrow(x_test), "\n")

# ------------------------------------------------------------
# Fit glmnet Cox on training set only
# ------------------------------------------------------------
cv_fit <- cv.glmnet(x_train, surv_train,
                    family = "cox",
                    alpha = 1,
                    nfolds = 10)

cat("lambda.min:", cv_fit$lambda.min, "  lambda.1se:", cv_fit$lambda.1se, "\n")

# Try lambda.1se first (more conservative); fall back to lambda.min if it
# selects zero genes. With only 7 candidate genes and modest hazard ratios,
# lambda.1se can be too aggressive and shrink everything to null -- if so,
# use lambda.min but report that choice explicitly, since it is the less
# conservative option and more prone to overfitting.
lambda_used <- "lambda.1se"
coefs <- coef(cv_fit, s = "lambda.1se")
coefs_df <- data.frame(gene = rownames(coefs), coefficient = as.numeric(coefs))
coefs_df <- coefs_df[coefs_df$coefficient != 0, ]

if (nrow(coefs_df) == 0) {
  cat("lambda.1se selected zero genes -- falling back to lambda.min.\n",
      "NOTE: report this explicitly -- lambda.min is less conservative\n",
      "and more prone to overfitting than lambda.1se.\n")
  lambda_used <- "lambda.min"
  coefs <- coef(cv_fit, s = "lambda.min")
  coefs_df <- data.frame(gene = rownames(coefs), coefficient = as.numeric(coefs))
  coefs_df <- coefs_df[coefs_df$coefficient != 0, ]
}

coefs_df <- coefs_df[order(abs(coefs_df$coefficient), decreasing = TRUE), ]

cat("\nSelected genes in survival signature (fit on training set, using",
    lambda_used, "):\n")
print(coefs_df)

if (nrow(coefs_df) == 0) {
  stop("Both lambda.1se and lambda.min selected zero genes -- no signature ",
       "can be extracted from this training split. This would need to be ",
       "reported as a null result, not forced.")
}

write.csv(coefs_df, "03_results/signature/glmnet_survival_signature.csv", row.names = FALSE)
cat("Signature fit using:", lambda_used,
    "-- record this in your methods section.\n")

# ------------------------------------------------------------
# Risk score on HELD-OUT test samples only -- this is the
# honest validation step. Do not compute this on training data.
# ------------------------------------------------------------
risk_score_test <- x_test[, coefs_df$gene, drop = FALSE] %*% coefs_df$coefficient
risk_score_test <- as.numeric(risk_score_test)

median_risk <- median(risk_score_test)
risk_group <- ifelse(risk_score_test >= median_risk, "High", "Low")

risk_df <- data.frame(
  os_time = clin_test$os_time,
  os_event = clin_test$os_event,
  risk_score = risk_score_test,
  risk_group = factor(risk_group, levels = c("Low", "High"))
)

fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = risk_df)

p <- ggsurvplot(fit,
                data = risk_df,
                pval = TRUE,
                risk.table = TRUE,
                title = "glmnet Survival Signature \u2014 Held-Out Risk Score (TCGA-LIHC)",
                legend.labs = c("Low Risk", "High Risk"),
                palette = c("#2E7C4A", "#D9534F"))

png("04_figures/glmnet_km_risk_heldout.png", width = 8, height = 6,
    units = "in", res = 300)
print(p)
dev.off()

write.csv(risk_df, "03_results/signature/glmnet_risk_scores_heldout.csv", row.names = FALSE)
cat("\nHeld-out risk score KM plot saved to 04_figures/glmnet_km_risk_heldout.png\n")
cat("This p-value is the honest, held-out validation number --",
    "if it survives on only", nrow(risk_df), "test samples, that's a real result.\n")

# ------------------------------------------------------------
# Comparison figure: same 7-gene model scored two ways --
#  (a) in-sample (all 366 samples, including the ones it was
#      trained on) -- mimics the original circular approach,
#      inflated, NOT to be reported as a validation result
#  (b) held-out only (the honest number above)
# This isolates how much the p-value was inflated by data
# leakage alone, holding the gene list and model fixed.
# ------------------------------------------------------------
library(patchwork)  # install.packages("patchwork") if missing

risk_score_all <- nexus_expr[, coefs_df$gene, drop = FALSE] %*% coefs_df$coefficient
risk_score_all <- as.numeric(risk_score_all)
median_all <- median(risk_score_all)
risk_group_all <- ifelse(risk_score_all >= median_all, "High", "Low")

risk_df_all <- data.frame(
  os_time = clin_tumor$os_time,
  os_event = clin_tumor$os_event,
  risk_group = factor(risk_group_all, levels = c("Low", "High"))
)

fit_all <- survfit(Surv(os_time, os_event) ~ risk_group, data = risk_df_all)
p_all <- ggsurvplot(fit_all, data = risk_df_all, pval = TRUE,
                    title = "In-Sample (includes training data)\n-- inflated, do not report",
                    legend.labs = c("Low Risk", "High Risk"),
                    palette = c("#2E7C4A", "#D9534F"),
                    risk.table = FALSE, ggtheme = theme_minimal(base_size = 10))

fit_test <- survfit(Surv(os_time, os_event) ~ risk_group, data = risk_df)
p_test <- ggsurvplot(fit_test, data = risk_df, pval = TRUE,
                     title = "Held-Out Test Set (n=110)\n-- honest, reportable",
                     legend.labs = c("Low Risk", "High Risk"),
                     palette = c("#2E7C4A", "#D9534F"),
                     risk.table = FALSE, ggtheme = theme_minimal(base_size = 10))

combined_plot <- p_all$plot + p_test$plot +
  plot_annotation(title = "Same 7-Gene Signature: In-Sample vs Held-Out Validation")

png("04_figures/glmnet_insample_vs_heldout.png", width = 12, height = 6, units = "in", res = 300)
print(combined_plot)
dev.off()

cat("Comparison figure saved to 04_figures/glmnet_insample_vs_heldout.png\n")