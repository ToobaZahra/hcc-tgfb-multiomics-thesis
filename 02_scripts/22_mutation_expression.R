# 22_mutation_expression.R
library(dplyr)

mut <- read.csv("03_results/mutation_nexus_genes.csv")
dge <- readRDS("03_results/dge_tcga_tumor_vs_normal.rds")

mut_summary <- mut %>%
  select(Hugo_Symbol, MutatedSamples) %>%
  arrange(desc(MutatedSamples))

dge_nexus <- dge[!is.na(dge$symbol) & dge$symbol %in% mut_summary$Hugo_Symbol,
                 c("symbol", "log2FoldChange", "padj")]

combined <- mut_summary %>%
  left_join(dge_nexus, by = c("Hugo_Symbol" = "symbol"))

write.csv(combined, "03_results/mutation_expression_cross.csv", row.names = FALSE)
cat("Mutation-expression cross:\n")
print(combined)