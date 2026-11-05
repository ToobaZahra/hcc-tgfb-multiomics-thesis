# sig_degs_export.R
# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/sig_degs_export.R
# Output : 03_results/dge/sig_DEGs_all_cohorts.xlsx (Summary + 12 DEG sheets)
# Cutoff : padj < 0.05 AND |log2FC| > LFC_CUT. Up = log2FC > 0, Down = log2FC < 0
# Dedup  : one row per gene symbol (best padj kept)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

dge_dir   <- "03_results/dge"
padj_cut  <- 0.05
LFC_CUT   <- 1        # change to 0.58 for 1.5-fold
COLLAPSE  <- TRUE     # FALSE = keep all probes/rows

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
  
  fc_col   <- pick(df, c("log2FoldChange", "logFC", "log2FC"), "log2FC", f)
  padj_col <- pick(df, c("padj", "adj.P.Val", "adj.P.value", "FDR"), "adjusted p-value", f)
  sym_col  <- pick(df, c("symbol", "Symbol", "SYMBOL", "Gene.symbol", "Gene.Symbol",
                         "gene_symbol", "GENE_SYMBOL", "Gene Symbol"), "gene symbol", f)
  
  df <- df[!is.na(df[[padj_col]]) & !is.na(df[[fc_col]]), ]
  n_tested_rows <- nrow(df)
  
  if (COLLAPSE) {
    s <- trimws(as.character(df[[sym_col]]))
    df <- df[!is.na(s) & s != "" & s != "NA", ]
    df <- df[order(df[[padj_col]]), ]
    df <- df[!duplicated(trimws(as.character(df[[sym_col]]))), ]
  }
  
  sig_padj <- df[df[[padj_col]] < padj_cut, ]
  n_padj   <- nrow(sig_padj)
  n_058    <- sum(abs(sig_padj[[fc_col]]) > 0.58)
  
  sig  <- sig_padj[abs(sig_padj[[fc_col]]) > LFC_CUT, ]
  sig  <- sig[order(sig[[padj_col]]), ]
  up   <- sig[sig[[fc_col]] > 0, ]
  down <- sig[sig[[fc_col]] < 0, ]
  up   <- up[order(-up[[fc_col]]), ]
  down <- down[order(down[[fc_col]]), ]
  
  for (s in list(list(paste0(co, "_all_sig"), sig),
                 list(paste0(co, "_up"),      up),
                 list(paste0(co, "_down"),    down))) {
    addWorksheet(wb, s[[1]])
    writeData(wb, s[[1]], s[[2]])
    freezePane(wb, s[[1]], firstRow = TRUE)
  }
  
  summary_df <- rbind(summary_df, data.frame(
    Cohort = co,
    Unique_genes_tested = nrow(df),
    Sig_padj_only = n_padj,
    Sig_padj_and_lfc0.58 = n_058,
    Sig_total = nrow(sig), Up = nrow(up), Down = nrow(down),
    Cutoff = paste0("padj < ", padj_cut, " & |log2FC| > ", LFC_CUT)))
  cat(sprintf("%-9s padj-only: %d | +|lfc|>0.58: %d | +|lfc|>%s: %d (up %d / down %d)\n",
              co, n_padj, n_058, LFC_CUT, nrow(sig), nrow(up), nrow(down)))
}

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
worksheetOrder(wb) <- c(length(names(wb)), seq_len(length(names(wb)) - 1))

out <- file.path(dge_dir, "sig_DEGs_all_cohorts.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
cat("Saved:", out, "\n")