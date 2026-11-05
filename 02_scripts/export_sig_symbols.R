# export_sig_full.R
# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/export_sig_full.R
# Reads the *_all_sig sheets from sig_DEGs_nominalP_all_cohorts.xlsx
# (nominal p < 0.05 & |log2FC| > 1, up + down combined, one row per gene)
# Writes one full-table CSV per cohort to 03_results/dge/network_inputs/

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
library(openxlsx)

xlsx <- "03_results/dge/sig_DEGs_nominalP_all_cohorts.xlsx"
outd <- "03_results/dge/network_inputs"
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

for (co in c("TCGA", "GSE14520", "GSE36376", "GSE76427")) {
  d <- read.xlsx(xlsx, sheet = paste0(co, "_all_sig"))
  fc <- intersect(c("log2FoldChange", "logFC", "log2FC"), names(d))[1]
  if (is.na(fc)) stop("No log2FC column in ", co, ". Columns: ", paste(names(d), collapse = ", "))
  d$direction <- ifelse(d[[fc]] > 0, "Up", "Down")
  
  out <- file.path(outd, sprintf("%s_sig_DEGs_full_nominalP.csv", co))
  write.csv(d, out, row.names = FALSE)
  cat(sprintf("%-9s %d genes (up %d / down %d) -> %s\n", co, nrow(d),
              sum(d$direction == "Up"), sum(d$direction == "Down"), out))
}


