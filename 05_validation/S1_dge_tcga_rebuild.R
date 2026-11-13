## S1_dge_tcga_rebuild.R
## Independent rebuild of TCGA-LIHC tumor vs normal DE, from the RAW GDC STAR count files.
## Save as: E:/MS-Project/hcc-tgfb-multiomics-thesis/05_validation/S1_dge_tcga_rebuild.R
## Run in the RStudio Console:  source("05_validation/S1_dge_tcga_rebuild.R")
## Outputs (05_validation/S1_dge_tcga/): rebuilt results, rebuilt VST, sample table, compare_report.txt, discordant genes.
## Original files are only READ, never modified.

suppressPackageStartupMessages({
  library(data.table); library(SummarizedExperiment); library(DESeq2); library(dplyr)
})
setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out <- "05_validation/S1_dge_tcga"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

LFC_TOL <- 0.01   # flag genes whose log2FC differs by more than this
log_lines <- character()
say <- function(...) { m <- paste0(...); cat(m, "\n"); log_lines <<- c(log_lines, m) }

## ============ S0: raw files vs saved SummarizedExperiment ============
rna_dir <- "data/raw/TCGA-LIHC/Transcriptome_Profiling/Gene_Expression_Quantification"
files <- list.files(rna_dir, pattern = "augmented_star_gene_counts\\.tsv$",
                    recursive = TRUE, full.names = TRUE)
say("== S0: raw RNA-seq files: ", length(files))

read_one <- function(f) {
  d <- fread(f, skip = "gene_id", select = c("gene_id", "gene_name", "unstranded"))
  d[!grepl("^N_", gene_id)]
}
first <- read_one(files[1])
raw_mat <- matrix(NA_real_, nrow(first), length(files),
                  dimnames = list(first$gene_id, basename(dirname(files))))
for (i in seq_along(files)) {
  d <- read_one(files[i])
  stopifnot(identical(d$gene_id, first$gene_id))
  raw_mat[, i] <- d$unstranded
}

se <- readRDS("data/raw/tcga_lihc_rnaseq_se.rds")
an <- if ("unstranded" %in% assayNames(se)) "unstranded" else assayNames(se)[1]
se_mat <- assay(se, an)
say("Saved SE: ", nrow(se_mat), " genes x ", ncol(se_mat), " samples (assay: ", an, ")")

common <- intersect(rownames(raw_mat), rownames(se_mat))
say("Genes raw / SE / common: ", nrow(raw_mat), " / ", nrow(se_mat), " / ", length(common))
rm_ <- raw_mat[common, , drop = FALSE]
sm_ <- se_mat[common, , drop = FALSE]
key <- function(m) sprintf("%.0f_%.0f", colSums(m), colSums(m * seq_len(nrow(m))))
kr <- key(rm_); ks <- key(sm_)
mt <- match(kr, ks)
ok <- !is.na(mt)
say("Raw files matched to an SE sample: ", sum(ok), " / ", length(kr))
say("SE samples with NO raw file:       ", sum(!(ks %in% kr)))
exact <- vapply(which(ok), function(i) all(rm_[, i] == sm_[, mt[i]]), logical(1))
say("Exact count-for-count equality:    ", sum(exact), " / ", sum(ok))
if (sum(ok) != length(kr) || !all(exact))
  say("  >>> WARNING: saved SE does not fully match raw files. Investigate before trusting any downstream result.")

## ============ S1: independent DESeq2 rebuild ============
meta <- data.frame(file_id = colnames(rm_)[ok],
                   barcode = colnames(se)[mt[ok]], stringsAsFactors = FALSE)
meta$patient <- substr(meta$barcode, 1, 12)
meta$type    <- substr(meta$barcode, 14, 15)
tt <- table(meta$type)
say("\n== S1: sample-type codes in the data: ", paste(names(tt), tt, sep = "=", collapse = ", "),
    "  (01 primary tumor, 02 recurrent, 11 solid normal)")

meta <- meta[meta$type %in% c("01", "11"), ]
meta <- meta[order(meta$barcode), ]
meta <- meta[!duplicated(paste(meta$patient, meta$type)), ]   # one vial per patient per type
meta$condition <- factor(ifelse(meta$type == "11", "normal", "tumor"), levels = c("normal", "tumor"))
say("Rebuild uses: ", sum(meta$condition == "tumor"), " tumor, ", sum(meta$condition == "normal"), " normal")
write.csv(meta, file.path(out, "rebuild_sample_table.csv"), row.names = FALSE)

counts <- round(rm_[, meta$file_id, drop = FALSE]); colnames(counts) <- meta$barcode
cd <- data.frame(condition = meta$condition, row.names = meta$barcode)
dds <- DESeqDataSetFromMatrix(counts, cd, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= min(table(cd$condition))
dds <- dds[keep, ]
say("Genes after filter (>=10 counts in >= ", min(table(cd$condition)), " samples): ", nrow(dds))
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "normal"), alpha = 0.05)
res_df <- data.frame(ensembl = rownames(res),
                     symbol = first$gene_name[match(rownames(res), first$gene_id)],
                     as.data.frame(res), stringsAsFactors = FALSE)
saveRDS(res_df, file.path(out, "dge_rebuilt.rds"))
write.csv(res_df, file.path(out, "dge_rebuilt.csv"), row.names = FALSE)
saveRDS(vst(dds), file.path(out, "vsd_rebuilt.rds"))     # reused by later stages
saveRDS(dds, file.path(out, "dds_rebuilt.rds"))

## ============ Compare with your original ============
say("\n== COMPARISON with original ==")
dds_o <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")
say("Original dds: ", ncol(dds_o), " samples, ", nrow(dds_o), " genes, design: ",
    paste(as.character(design(dds_o)), collapse = " "))
say("Original sample overlap: in both = ", length(intersect(colnames(dds_o), colnames(dds))),
    "; only original = ", length(setdiff(colnames(dds_o), colnames(dds))),
    "; only rebuild = ", length(setdiff(colnames(dds), colnames(dds_o))))
if (length(setdiff(colnames(dds_o), colnames(dds))) > 0) {
  o_only <- setdiff(colnames(dds_o), colnames(dds))
  say("  Sample types of original-only samples: ",
      paste(names(table(substr(o_only, 14, 15))), table(substr(o_only, 14, 15)), sep = "=", collapse = ", "))
}
say("Gene overlap (filtered sets): both = ", length(intersect(rownames(dds_o), rownames(dds))),
    "; only original = ", length(setdiff(rownames(dds_o), rownames(dds))),
    "; only rebuild = ", length(setdiff(rownames(dds), rownames(dds_o))))

op <- "03_results/dge/dge_tcga_tumor_vs_normal.rds"
if (!file.exists(op)) op <- "03_results/dge_tcga_tumor_vs_normal.rds"
orig <- as.data.frame(readRDS(op))
say("Original results file: ", op)
say("Original columns: ", paste(colnames(orig), collapse = ", "))

pick <- function(nms, pats) { for (p in pats) { h <- grep(p, nms, ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }; NA }
c_sym <- pick(colnames(orig), c("^symbol$", "gene_symbol", "^gene_name$", "^gene$", "SYMBOL"))
c_lfc <- pick(colnames(orig), c("log2FoldChange", "log2FC", "logFC"))
c_pad <- pick(colnames(orig), c("^padj$", "adj\\.P\\.Val", "FDR"))
say("Matching on: symbol = ", c_sym, ", lfc = ", c_lfc, ", padj = ", c_pad)
if (is.na(c_sym)) { orig$symbol_tmp <- rownames(orig); c_sym <- "symbol_tmp" }

o <- data.frame(symbol = orig[[c_sym]], lfc_o = orig[[c_lfc]], padj_o = orig[[c_pad]], stringsAsFactors = FALSE)
r <- data.frame(symbol = res_df$symbol, lfc_r = res_df$log2FoldChange, padj_r = res_df$padj, stringsAsFactors = FALSE)
o <- o[!is.na(o$symbol) & o$symbol != "", ]; r <- r[!is.na(r$symbol) & r$symbol != "", ]
o <- o[order(o$padj_o), ]; o <- o[!duplicated(o$symbol), ]    # one row per symbol (best padj)
r <- r[order(r$padj_r), ]; r <- r[!duplicated(r$symbol), ]
mm <- inner_join(o, r, by = "symbol")
mm$dlfc <- mm$lfc_r - mm$lfc_o
say("\nSymbols: original ", nrow(o), ", rebuild ", nrow(r), ", common ", nrow(mm))
say("log2FC: Pearson r = ", round(cor(mm$lfc_o, mm$lfc_r, use = "complete.obs"), 5),
    "; max |diff| = ", round(max(abs(mm$dlfc), na.rm = TRUE), 4),
    "; genes with |diff| > ", LFC_TOL, " = ", sum(abs(mm$dlfc) > LFC_TOL, na.rm = TRUE))

sigA_o <- mm$padj_o < 0.05;  sigA_r <- mm$padj_r < 0.05
sigB_o <- sigA_o & abs(mm$lfc_o) > 1; sigB_r <- sigA_r & abs(mm$lfc_r) > 1
sigC_o <- mm$padj_o < 0.1;   sigC_r <- mm$padj_r < 0.1
cnt <- function(s, l) sprintf("up=%d down=%d", sum(s & l > 0, na.rm = TRUE), sum(s & l < 0, na.rm = TRUE))
say("padj<0.05           original ", cnt(sigA_o, mm$lfc_o), " | rebuild ", cnt(sigA_r, mm$lfc_r))
say("padj<0.05 & |lfc|>1 original ", cnt(sigB_o, mm$lfc_o), " | rebuild ", cnt(sigB_r, mm$lfc_r))
say("padj<0.1  (DESeq2 default; your reported 15,991 up / 6,058 down) original ",
    cnt(sigC_o, mm$lfc_o), " | rebuild ", cnt(sigC_r, mm$lfc_r))
say("Significance calls that differ: padj<0.05 -> ", sum(sigA_o != sigA_r, na.rm = TRUE),
    "; padj<0.05&|lfc|>1 -> ", sum(sigB_o != sigB_r, na.rm = TRUE))
say("Direction of effect flipped (sig in both): ", sum(sigA_o & sigA_r & sign(mm$lfc_o) != sign(mm$lfc_r), na.rm = TRUE))

disc <- mm[abs(mm$dlfc) > LFC_TOL | (sigB_o != sigB_r), ]
write.csv(disc[order(-abs(disc$dlfc)), ], file.path(out, "discordant_genes.csv"), row.names = FALSE)
say("Discordant genes written: ", nrow(disc), " (discordant_genes.csv)")

## ---- the 26 nexus genes ----
nf <- c("03_results/dge/dge_nexus_genes.csv", "03_results/dge_nexus_genes.csv")
nf <- nf[file.exists(nf)][1]
if (!is.na(nf)) {
  nx <- read.csv(nf, stringsAsFactors = FALSE)
  ns <- pick(colnames(nx), c("^symbol$", "gene_symbol", "^gene_name$", "^gene$", "SYMBOL"))
  syms <- if (!is.na(ns)) unique(nx[[ns]]) else NULL
  if (!is.null(syms)) {
    tab <- mm[mm$symbol %in% syms, ]
    tab$sig_o <- tab$padj_o < 0.05 & abs(tab$lfc_o) > 1
    tab$sig_r <- tab$padj_r < 0.05 & abs(tab$lfc_r) > 1
    write.csv(tab, file.path(out, "nexus_genes_compare.csv"), row.names = FALSE)
    say("\nNexus genes found in both: ", nrow(tab), " of ", length(syms),
        "; max |dlfc| = ", round(max(abs(tab$dlfc)), 4),
        "; sig call differs for ", sum(tab$sig_o != tab$sig_r), " gene(s)")
    say("Missing from rebuild: ", paste(setdiff(syms, mm$symbol), collapse = ", "))
  }
}

writeLines(log_lines, file.path(out, "compare_report.txt"))
cat("\nDONE. Paste the full console output (or compare_report.txt) back to Claude.\n")