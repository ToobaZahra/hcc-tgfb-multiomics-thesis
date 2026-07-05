# 28_survival.R
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)

# Load data
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds")
res_df <- readRDS("03_results/dge/dge_tcga_tumor_vs_normal.rds")

# Build survival object
clin$os_event <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$os_time <- ifelse(!is.na(clin$days_to_death),
                       clin$days_to_death,
                       clin$days_to_last_follow_up)

clin_surv <- clin[!is.na(clin$os_time) & clin$os_time > 0, ]

# Match with expression
vsd$patient <- substr(colnames(vsd), 1, 12)
common <- intersect(clin_surv$submitter_id, vsd$patient)
clin_matched <- clin_surv[clin_surv$submitter_id %in% common, ]
vsd_matched <- vsd[, vsd$patient %in% common]
tumor_keep <- vsd_matched$sample_type == "Primary Tumor"
vsd_tumor <- vsd_matched[, tumor_keep]
clin_tumor <- clin_matched[match(vsd_tumor$patient, clin_matched$submitter_id), ]

cat("Samples:", ncol(vsd_tumor), "\n")
cat("Events:", sum(clin_tumor$os_event), "\n")

# Univariate Cox for all nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

expr_mat <- assay(vsd_tumor)
cox_results <- data.frame()

for (i in 1:nrow(gene_map)) {
  sym <- gene_map$symbol[i]
  ens <- gene_map$ensembl_id[i]
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) == 0) next
  expr <- expr_mat[row_match[1], ]
  df <- data.frame(os_time = clin_tumor$os_time,
                   os_event = clin_tumor$os_event,
                   expr = as.numeric(expr))
  df <- df[!is.na(df$os_time), ]
  fit <- coxph(Surv(os_time, os_event) ~ expr, data = df)
  s <- summary(fit)
  cox_results <- rbind(cox_results, data.frame(
    gene = sym,
    HR = round(s$coefficients[1, "exp(coef)"], 3),
    CI_low = round(s$conf.int[1, "lower .95"], 3),
    CI_high = round(s$conf.int[1, "upper .95"], 3),
    pval = round(s$coefficients[1, "Pr(>|z|)"], 4)
  ))
}

cox_results$padj <- p.adjust(cox_results$pval, method = "BH")
cox_results <- cox_results[order(cox_results$pval), ]

write.csv(cox_results, "03_results/survival/cox_univariate.csv", row.names = FALSE)
cat("Cox results:\n")
print(cox_results)