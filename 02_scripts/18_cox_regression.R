# 18_cox_regression.R
library(survival)
library(SummarizedExperiment)
library(dplyr)

# Load data
clin_tumor <- readRDS("03_results/survival_clinical_matched.rds")
vsd_tumor <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds")
res_df <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

expr_mat <- assay(vsd_tumor)

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR2","SMAD3","SMAD4","MYC","KLF4",
                 "DUSP1","HDAC11","SNAI1","IL6")

# Match gene symbols to Ensembl IDs
gene_map <- res_df[!is.na(res_df$symbol) & res_df$symbol %in% nexus_genes,
                   c("symbol", "ensembl_id")]
gene_map <- gene_map[!duplicated(gene_map$symbol), ]

# Run univariate Cox for each gene
cox_results <- data.frame()

for (i in 1:nrow(gene_map)) {
  sym <- gene_map$symbol[i]
  ens <- gene_map$ensembl_id[i]
  
  row_match <- grep(paste0("^", ens), rownames(expr_mat))
  if (length(row_match) == 0) next
  
  expr <- expr_mat[row_match[1], ]
  
  df <- data.frame(
    os_time  = clin_tumor$os_time,
    os_event = clin_tumor$os_event,
    expr     = as.numeric(expr)
  )
  df <- df[!is.na(df$os_time), ]
  
  fit <- coxph(Surv(os_time, os_event) ~ expr, data = df)
  s <- summary(fit)
  
  cox_results <- rbind(cox_results, data.frame(
    gene   = sym,
    HR     = round(s$coefficients[1, "exp(coef)"], 3),
    CI_low = round(s$conf.int[1, "lower .95"], 3),
    CI_high= round(s$conf.int[1, "upper .95"], 3),
    pval   = round(s$coefficients[1, "Pr(>|z|)"], 4)
  ))
}

cox_results <- cox_results[order(cox_results$pval), ]
write.csv(cox_results, "03_results/cox_nexus_genes.csv", row.names = FALSE)

cat("Cox regression results:\n")
print(cox_results)