# 32_signature_evaluation.R
library(glmnet)
library(pROC)

cv_fit <- readRDS("03_results/signature/glmnet_fit.rds")

# Predict on same data (in-sample evaluation)
pred_prob <- predict(cv_fit, newx = nexus_expr, 
                     s = "lambda.min", type = "response")

# ROC curve
roc_obj <- roc(as.numeric(y == "Late"), as.numeric(pred_prob))
auc_val <- auc(roc_obj)
cat("AUC:", round(auc_val, 3), "\n")

# Plot ROC
png("04_figures/roc_stage_signature.png", width = 6, height = 6, 
    units = "in", res = 300)
plot(roc_obj, col = "#B2182B", lwd = 2,
     main = paste0("ROC — Stage Signature (AUC = ", round(auc_val, 3), ")"))
dev.off()
cat("ROC plot saved.\n")