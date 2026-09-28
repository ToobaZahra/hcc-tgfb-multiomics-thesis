# 21_mutation_analysis.R

library(TCGAbiolinks)
library(maftools)

# ============================================================
# 1. Query TCGA-LIHC somatic mutation data
# ============================================================

query_maf <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(
  query_maf,
  method = "api",
  directory = "data/raw/GDCdata"
)

maf_data <- GDCprepare(
  query_maf,
  directory = "data/raw/GDCdata"
)

saveRDS(
  maf_data,
  "data/raw/tcga_lihc_maf.rds"
)

cat("MAF data saved.\n")
cat("Dimensions:", dim(maf_data), "\n")


# ============================================================
# 2. Convert data to maftools MAF object
# ============================================================

maf <- read.maf(maf = maf_data)

cat("\nMAF Summary:\n")
print(getSampleSummary(maf)[1:10, ])


# ============================================================
# 3. Define TGF-β nexus genes
# ============================================================

nexus_genes <- c(
  "TGFB1", "TGFBR1", "TGFBR2",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  "MAPK8", "MAPK9", "MAPK10",
  "DUSP1", "DUSP4", "DUSP10",
  "MYC", "CDKN1A", "CDKN2B",
  "SNAI1", "TNF", "IL6", "IL10", "IL37",
  "HDAC11", "NPC1", "CCDC110",
  "TGFBRAP1", "KLF4"
)


# ============================================================
# 4. Generate oncoplot for nexus genes
# ============================================================

png(
  "04_figures/oncoplot_nexus.png",
  width = 14,
  height = 8,
  units = "in",
  res = 300
)

oncoplot(
  maf = maf,
  genes = nexus_genes,
  titleText = "TGF-β Nexus Genes — Somatic Mutations (TCGA-LIHC)"
)

dev.off()

cat("Oncoplot saved.\n")


# ============================================================
# 5. Calculate mutation frequency of nexus genes
# ============================================================

mut_summary <- getGeneSummary(maf)

nexus_mut_summary <- mut_summary[
  mut_summary$Hugo_Symbol %in% nexus_genes,
]

nexus_mut_summary <- nexus_mut_summary[
  order(
    nexus_mut_summary$MutatedSamples,
    decreasing = TRUE
  ),
]


# ============================================================
# 6. Save mutation summary
# ============================================================

write.csv(
  nexus_mut_summary,
  "03_results/mutation_nexus_genes.csv",
  row.names = FALSE
)

cat("\nNexus gene mutation frequencies:\n")

print(
  nexus_mut_summary[
    ,
    c("Hugo_Symbol", "MutatedSamples", "AlteredSamples")
  ]
)

cat("\nMutation analysis completed successfully.\n")