# 21_mutation_analysis.R
library(TCGAbiolinks)
library(maftools)

# Query somatic mutations
query_maf <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query_maf, method = "api", 
            directory = "data/raw/GDCdata")

maf_data <- GDCprepare(query_maf, 
                       directory = "data/raw/GDCdata")

saveRDS(maf_data, "data/raw/tcga_lihc_maf.rds")
cat("MAF data saved.\n")
cat("Dimensions:", dim(maf_data), "\n")



# Load maftools
library(maftools)

# Convert to MAF object
maf <- read.maf(maf = maf_data)

# Summary
cat("\nMAF Summary:\n")
print(getSampleSummary(maf)[1:10, ])

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# Oncoprint for nexus genes
png("04_figures/oncoplot_nexus.png", width = 14, height = 8, 
    units = "in", res = 300)
oncoplot(maf = maf, genes = nexus_genes, 
         title = "TGF-β Nexus Genes — Somatic Mutations (TCGA-LIHC)")
dev.off()
cat("Oncoplot saved.\n")

# Mutation frequency in nexus genes
nexus_mut <- subsetMaf(maf, genes = nexus_genes)
mut_summary <- getGeneSummary(maf)
nexus_mut_summary <- mut_summary[mut_summary$Hugo_Symbol %in% nexus_genes, ]
nexus_mut_summary <- nexus_mut_summary[order(nexus_mut_summary$MutatedSamples, 
                                             decreasing = TRUE), ]

write.csv(nexus_mut_summary, "03_results/mutation_nexus_genes.csv", 
          row.names = FALSE)

cat("\nNexus gene mutation frequencies:\n")
print(nexus_mut_summary[, c("Hugo_Symbol", "MutatedSamples", 
                            "AlteredSamples")])

dev.off()  # close the png device that opened

png("04_figures/oncoplot_nexus.png", width = 14, height = 8, 
    units = "in", res = 300)
oncoplot(maf = maf, genes = nexus_genes)
dev.off()
cat("Oncoplot saved.\n")

mut_summary <- getGeneSummary(maf)
nexus_mut_summary <- mut_summary[mut_summary$Hugo_Symbol %in% nexus_genes, ]
nexus_mut_summary <- nexus_mut_summary[order(nexus_mut_summary$MutatedSamples, 
                                             decreasing = TRUE), ]

write.csv(nexus_mut_summary, "03_results/mutation_nexus_genes.csv", 
          row.names = FALSE)

cat("\nNexus gene mutation frequencies:\n")
print(nexus_mut_summary[, c("Hugo_Symbol", "MutatedSamples", "AlteredSamples")])