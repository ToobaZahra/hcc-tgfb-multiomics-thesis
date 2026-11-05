## 63: Merge the cross-cohort tables of DUSP10 (file 23) and SMAD2, MYC, MAPK8, HDAC11 (files 32_*) into one file.
## Run AFTER scripts 53 (DUSP10) and 62 (other genes). Output: one CSV (and an Excel file with one sheet per gene, if openxlsx is installed).
suppressPackageStartupMessages(library(dplyr))

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"

rd <- function(f, gene) {
  p <- file.path(a, f)
  if (!file.exists(p)) { cat("Missing:", f, "\n"); return(NULL) }
  x <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE); x$gene <- gene; x }

parts <- list(
  rd("23_DUSP10_cross_cohort_summary.csv",  "DUSP10"),
  rd("32_SMAD2_cross_cohort_summary.csv",   "SMAD2"),
  rd("32_MYC_cross_cohort_summary.csv",     "MYC"),
  rd("32_MAPK8_cross_cohort_summary.csv",   "MAPK8"),
  rd("32_HDAC11_cross_cohort_summary.csv",  "HDAC11"))

cohorts <- c("TCGA-LIHC", "GSE14520", "GSE76427", "GSE56588", "GSE54503")
merged <- bind_rows(parts)
for (cc in setdiff(cohorts, names(merged))) merged[[cc]] <- NA_character_
merged <- merged %>%
  mutate(gene = factor(gene, c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")),
         p_type = ifelse(grepl("padj", analysis), "adjusted (padj)", "raw p"),
         probe_rows = grepl("^Methylation, probe", layer)) %>%
  arrange(gene, probe_rows, layer, analysis) %>%
  dplyr::select(gene, layer, analysis, statistic, all_of(cohorts), p_type)

write.csv(merged, file.path(a, "33_all_genes_cross_cohort_merged.csv"), row.names = FALSE, na = "NA")
cat("Merged table:", nrow(merged), "rows for genes:", paste(unique(as.character(merged$gene)), collapse = ", "), "\n")
print(table(merged$gene))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All genes"); openxlsx::writeData(wb, "All genes", merged)
  for (gn in levels(merged$gene)) {
    s <- merged[merged$gene == gn, ]
    if (nrow(s) == 0) next
    openxlsx::addWorksheet(wb, gn); openxlsx::writeData(wb, gn, s) }
  note <- data.frame(Note = c(
    "p-values are raw unless the analysis name says padj.",
    "TCGA expression is in VST units and GEO expression is log2: compare direction, not size.",
    "NA = analysis not possible in that cohort (expression cohorts have no methylation; GSE56588 and GSE54503 have methylation only).",
    "Survival: HR per 1 SD. TCGA from the patient-level clinical table.",
    "DUSP10 probe rows are cg23932332; for the other genes the probe is the one with the largest |rho| with expression (selected on TCGA, so optimistic).",
    "GSE54503: paired tumor vs adjacent tissue (60 pairs)."))
  openxlsx::addWorksheet(wb, "Notes"); openxlsx::writeData(wb, "Notes", note)
  openxlsx::saveWorkbook(wb, file.path(a, "33_all_genes_cross_cohort_merged.xlsx"), overwrite = TRUE)
  cat("Saved 33_all_genes_cross_cohort_merged.xlsx (sheets: All genes, one per gene, Notes)\n")
} else cat("Install openxlsx (install.packages('openxlsx')) if you also want the Excel version.\n")

cat("\nSaved 33_all_genes_cross_cohort_merged.csv\n")
cat("Notes: p-values are raw unless marked padj (see p_type). Compare direction across cohorts, not size.\n")