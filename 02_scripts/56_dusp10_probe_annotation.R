## 56: Genomic context of the DUSP10 promoter probes (esp. cg23932332)
## CpG island relation, gene region, regulatory feature, enhancer/DHS, from the Illumina 450K annotation.
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
key <- "cg23932332"

# ---- 1. Annotation: use the TCGA object if it carries it, else the Bioconductor 450K annotation ----
pt <- read.csv(file.path(a, "19_DUSP10_probe_table.csv")) %>% select(probe, dist_TSS, mean_beta, rho_expr, padj_expr)

pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (!requireNamespace(pkg, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c(pkg, "minfi"), update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(minfi); library(pkg, character.only = TRUE) })
ann <- as.data.frame(getAnnotation(get(pkg)))
want <- c("chr", "pos", "strand", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
          "Regulatory_Feature_Group", "Enhancer", "DHS", "Phantom5_Enhancers")
want <- intersect(want, colnames(ann))
an <- ann[pt$probe, want, drop = FALSE]; an$probe <- rownames(an)

res <- pt %>% left_join(an, by = "probe") %>% arrange(dist_TSS)
write.csv(res, file.path(a, "26_DUSP10_probe_annotation.csv"), row.names = FALSE)

cat("== DUSP10 promoter probes: genomic context ==\n")
print(res %>% select(probe, dist_TSS, any_of(c("Relation_to_Island", "UCSC_RefGene_Group",
                                               "Regulatory_Feature_Group", "Enhancer", "DHS"))), right = FALSE)
cat("\n== Key probe", key, "==\n"); print(t(res[res$probe == key, ]))

# ---- 2. Quick context: is the key probe different from the other 16? ----
cat("\nRelation to CpG island across the 17 probes:\n"); print(table(res$Relation_to_Island, useNA = "ifany"))
cat("\nNote: annotation is hg19 (450K manifest ilmn12). Island relations: Island, N_Shore, S_Shore, N_Shelf, S_Shelf, OpenSea.\n")
cat("Check that the strand and UCSC_RefGene_Group entries refer to DUSP10 (a probe can map to several genes).\n")