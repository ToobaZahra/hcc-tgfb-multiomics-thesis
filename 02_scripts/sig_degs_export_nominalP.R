# sig_degs_export_nominalP.R
# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/sig_degs_export_nominalP.R
# Output : 03_results/dge/sig_DEGs_nominalP_all_cohorts.xlsx (Summary + 12 DEG sheets)
# Cutoff : nominal p < P_CUT AND |log2FC| > LFC_CUT (NO multiple-testing correction)
# Dedup  : one row per gene symbol (lowest p kept)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

dge_dir  <- "03_results/dge"
P_CUT    <- 0.05
LFC_CUT  <- 0.5
COLLAPSE <- TRUE

cohorts <- list(
  TCGA     = "dge_tcga_tumor_vs_normal.csv",
  GSE14520 = "limma_GSE14520_tumor_vs_normal.csv",
  GSE36376 = "limma_GSE36376_tumor_vs_normal.csv",
  GSE76427 = "limma_GSE76427_tumor_vs_normal.csv"
)

pick <- function(df, candidates, what, f) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0)
    stop("No ", what, " column in ", f, ". Columns found: ",
         paste(names(df), collapse = ", "))
  hit[1]
}

wb <- createWorkbook()
summary_df <- data.frame()

for (co in names(cohorts)) {
  f  <- file.path(dge_dir, cohorts[[co]])
  df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  
  fc_col  <- pick(df, c("log2FoldChange", "logFC", "log2FC"), "log2FC", f)
  p_col   <- pick(df, c("pvalue", "P.Value", "PValue", "p.value"), "nominal p-value", f)
  sym_col <- pick(df, c("symbol", "Symbol", "SYMBOL", "Gene.symbol", "Gene.Symbol",
                        "gene_symbol", "GENE_SYMBOL", "Gene Symbol"), "gene symbol", f)
  
  df <- df[!is.na(df[[p_col]]) & !is.na(df[[fc_col]]), ]
  
  if (COLLAPSE) {
    s <- trimws(as.character(df[[sym_col]]))
    df <- df[!is.na(s) & s != "" & s != "NA", ]
    df <- df[order(df[[p_col]]), ]
    df <- df[!duplicated(trimws(as.character(df[[sym_col]]))), ]
  }
  
  sig_p <- df[df[[p_col]] < P_CUT, ]
  sig   <- sig_p[abs(sig_p[[fc_col]]) > LFC_CUT, ]
  sig   <- sig[order(sig[[p_col]]), ]
  up    <- sig[sig[[fc_col]] > 0, ]
  down  <- sig[sig[[fc_col]] < 0, ]
  up    <- up[order(-up[[fc_col]]), ]
  down  <- down[order(down[[fc_col]]), ]
  
  for (s in list(list(paste0(co, "_all_sig"), sig),
                 list(paste0(co, "_up"),      up),
                 list(paste0(co, "_down"),    down))) {
    addWorksheet(wb, s[[1]])
    writeData(wb, s[[1]], s[[2]])
    freezePane(wb, s[[1]], firstRow = TRUE)
  }
  
  summary_df <- rbind(summary_df, data.frame(
    Cohort = co, Unique_genes_tested = nrow(df),
    Sig_p_only = nrow(sig_p),
    Sig_total = nrow(sig), Up = nrow(up), Down = nrow(down),
    Cutoff = paste0("nominal p < ", P_CUT, " & |log2FC| > ", LFC_CUT)))
  cat(sprintf("%-9s p-only: %d | +|lfc|>%s: %d (up %d / down %d)\n",
              co, nrow(sig_p), LFC_CUT, nrow(sig), nrow(up), nrow(down)))
}

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
worksheetOrder(wb) <- c(length(names(wb)), seq_len(length(names(wb)) - 1))

out <- file.path(dge_dir, "sig_DEGs_nominalP_all_cohorts.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
cat("Saved:", out, "\n")