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

# ------------------------------------------------------------
# Nexus genes significant in univariate Cox (from 28_survival.R)
# Pulled programmatically instead of hardcoded, so this stays in
# sync if the univariate results ever change.
# ------------------------------------------------------------
cox_univariate <- read.csv("03_results/survival/cox_univariate.csv")
sig_genes <- cox_univariate$gene[cox_univariate$pval < 0.05]

cat("Genes selected for stage-adjusted Cox (pval < 0.05 in univariate):\n")
print(sig_genes)

# Carry forward each gene's univariate FDR status so it's visible in the
# stage-adjusted table too. Only 2/26 genes (DUSP10, NPC1) survive BH
# correction across all 26 nexus genes -- the rest are nominal (p<0.05,
# padj>0.05) hits. Report both numbers, don't call the nominal ones
# "significant" without qualifying which threshold you mean.
univariate_fdr_status <- setNames(
  ifelse(cox_univariate$padj < 0.05, "FDR-significant", "Nominal only"),
  cox_univariate$gene
)

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
cox_stage_results$univariate_fdr <- univariate_fdr_status[cox_stage_results$gene]
cox_stage_results <- cox_stage_results[order(cox_stage_results$pval), ]

write.csv(cox_stage_results, "03_results/survival/cox_stage_adjusted.csv", row.names = FALSE)
cat("Stage-adjusted Cox results:\n")
print(cox_stage_results)

cat("\nNote: 'univariate_fdr' reflects each gene's FDR status in the ORIGINAL",
    "26-gene univariate screen (28_survival.R), not in this stage-adjusted model.",
    "\nOnly", sum(univariate_fdr_status == "FDR-significant"),
    "/ 26 genes (DUSP10, NPC1) were FDR-significant in the univariate screen --",
    "the rest reported here were nominal (p<0.05) hits only.\n")

library(ggplot2)

# Merge univariate and stage-adjusted results for the 7 genes
univ_sub <- cox_univariate[cox_univariate$gene %in% cox_stage_results$gene,
                           c("gene", "HR", "CI_low", "CI_high")]
univ_sub$model <- "Univariate"

stage_sub <- cox_stage_results[, c("gene", "HR", "CI_low", "CI_high")]
stage_sub$model <- "Stage-adjusted"

combined <- rbind(univ_sub, stage_sub)

# order genes by stage-adjusted p-value (same order as your printed table)
gene_order <- cox_stage_results$gene[order(cox_stage_results$pval)]
combined$gene <- factor(combined$gene, levels = rev(gene_order))
combined$model <- factor(combined$model, levels = c("Univariate", "Stage-adjusted"))

p <- ggplot(combined, aes(x = HR, y = gene, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Univariate" = "#999999", "Stage-adjusted" = "#C0392B")) +
  labs(title = "Univariate vs Stage-Adjusted Cox Hazard Ratios",
       subtitle = "Nexus genes with nominal p<0.05 in univariate screen",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave("04_figures/cox_univariate_vs_stage_adjusted.png", p, width = 8, height = 6, dpi = 300)
print(p)