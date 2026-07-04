# 25_integrated_summary.R
library(dplyr)

# Load all results
dge <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")
cox <- read.csv("03_results/cox_nexus_genes.csv")
meth <- read.csv("03_results/dmeth_nexus_summary.csv")
muts <- read.csv("03_results/mutation_nexus_genes.csv")
variants <- read.csv("03_results/nexus_variants_summary.csv")
replication <- read.csv("03_results/nexus_replication_summary.csv")

nexus <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# DGE
dge_clean <- dge[dge$symbol %in% nexus & !is.na(dge$symbol), ] %>%
  select(symbol, logFC_TCGA = log2FoldChange, padj_TCGA = padj)

# Methylation
meth_clean <- meth %>%
  select(symbol = gene, n_meth_sig = n_sig, mean_delta_beta = mean_logFC)

# Cox
cox_clean <- cox %>%
  select(symbol = gene, HR, cox_pval = pval)

# Mutations
muts_clean <- muts %>%
  select(symbol = Hugo_Symbol, n_somatic_muts = MutatedSamples)

# Variants
var_clean <- variants %>%
  select(symbol = gene_symbol, n_pathogenic = n_pathogenic,
         n_high_impact = total_high_impact)

# Replication
rep_clean <- replication %>%
  select(symbol, replicated = consistent)

# Merge all
summary_table <- data.frame(symbol = nexus) %>%
  left_join(dge_clean, by = "symbol") %>%
  left_join(meth_clean, by = "symbol") %>%
  left_join(cox_clean, by = "symbol") %>%
  left_join(muts_clean, by = "symbol") %>%
  left_join(var_clean, by = "symbol") %>%
  left_join(rep_clean, by = "symbol")

write.csv(summary_table, "03_results/nexus_integrated_summary.csv", row.names = FALSE)
cat("Integrated summary table saved.\n")
print(summary_table)