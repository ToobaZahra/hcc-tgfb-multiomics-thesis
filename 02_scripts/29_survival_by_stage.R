# 29_survival_by_stage.R
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)

clin_tumor <- readRDS("03_results/survival/survival_clinical_matched.rds")
vsd_tumor <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds")
res_df <- readRDS("03_results/dge/dge_tcga_tumor_vs_normal.rds")

# Add stage group
clin_tumor$stage_group <- ifelse(
  clin_tumor$ajcc_pathologic_stage %in% c("Stage I", "Stage II"), "Early",
  ifelse(clin_tumor$ajcc_pathologic_stage %in% 
           c("Stage III","Stage IIIA","Stage IIIB","Stage IV"), "Late", NA))

cat("Stage distribution:\n")
print(table(clin_tumor$stage_group, useNA = "always"))



# Nexus genes significant in Cox
sig_genes <- c("DUSP10", "NPC1", "SMAD2", "HDAC11", "CDKN2B", "CCDC110")

gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% sig_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

expr_mat <- assay(vsd_tumor)
cox_stage_results <- data.frame()

for (i in 1:nrow(gene_map)) {
  sym <- gene_map$symbol[i]
  ens <- gene_map$ensembl_id[i]
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) == 0) next
  expr <- expr_mat[row_match[1], ]
  
  df <- data.frame(
    os_time = clin_tumor$os_time,
    os_event = clin_tumor$os_event,
    expr = as.numeric(expr),
    stage = clin_tumor$stage_group
  )
  df <- df[!is.na(df$os_time) & !is.na(df$stage), ]
  
  fit <- coxph(Surv(os_time, os_event) ~ expr + stage, data = df)
  s <- summary(fit)
  
  cox_stage_results <- rbind(cox_stage_results, data.frame(
    gene = sym,
    HR = round(s$coefficients[1, "exp(coef)"], 3),
    CI_low = round(s$conf.int[1, "lower .95"], 3),
    CI_high = round(s$conf.int[1, "upper .95"], 3),
    pval = round(s$coefficients[1, "Pr(>|z|)"], 4)
  ))
}

cox_stage_results$padj <- p.adjust(cox_stage_results$pval, method = "BH")
cox_stage_results <- cox_stage_results[order(cox_stage_results$pval), ]

write.csv(cox_stage_results, "03_results/survival/cox_stage_adjusted.csv", row.names = FALSE)
cat("Stage-adjusted Cox results:\n")
print(cox_stage_results)