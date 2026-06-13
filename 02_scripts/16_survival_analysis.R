# 16_survival_analysis.R
library(survival)
library(survminer)
library(DESeq2)
library(dplyr)

# Load clinical data
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")
cat("Clinical data dimensions:", dim(clin), "\n")
cat("Clinical columns:\n")
print(colnames(clin))

# Build survival object
# OS = overall survival
# Event: vital_status == "Dead" = 1, "Alive" = 0
clin$os_event <- ifelse(clin$vital_status == "Dead", 1, 0)

# OS time: days_to_death for dead, days_to_last_follow_up for alive
clin$os_time <- ifelse(!is.na(clin$days_to_death), 
                       clin$days_to_death, 
                       clin$days_to_last_follow_up)

# Remove missing
clin_surv <- clin[!is.na(clin$os_time) & clin$os_time > 0, ]
cat("Samples with survival data:", nrow(clin_surv), "\n")
cat("Events (deaths):", sum(clin_surv$os_event), "\n")

# Load VST expression
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds")

# Match samples
vsd$patient <- substr(colnames(vsd), 1, 12)
common <- intersect(clin_surv$submitter_id, vsd$patient)
cat("Matched samples:", length(common), "\n")

# Subset
clin_matched <- clin_surv[clin_surv$submitter_id %in% common, ]
vsd_matched <- vsd[, vsd$patient %in% common]

# Tumor only
tumor_keep <- vsd_matched$sample_type == "Primary Tumor"
vsd_tumor <- vsd_matched[, tumor_keep]
clin_tumor <- clin_matched[match(vsd_tumor$patient, clin_matched$submitter_id), ]

cat("Tumor samples for survival:", ncol(vsd_tumor), "\n")

# Save matched data
saveRDS(clin_tumor, "03_results/survival_clinical_matched.rds")
saveRDS(vsd_tumor, "data/raw/tcga_lihc_vsd_tumor_matched.rds")
cat("Survival data saved.\n")