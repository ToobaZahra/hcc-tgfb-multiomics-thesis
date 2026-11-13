## S0_probe.R -- READ-ONLY inspection of cached objects. Changes nothing except writing one text file.
## Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/05_validation/S0_probe.R
## Run in the RStudio Console:  source("05_validation/S0_probe.R")
## Then paste the contents of 05_validation/S0_probe_output.txt back into the chat.

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
dir.create("05_validation", showWarnings = FALSE)
sink("05_validation/S0_probe_output.txt", split = TRUE)

suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2) })
safe <- function(label, expr) {
  cat("\n====", label, "====\n")
  tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

safe("Session", {
  print(R.version.string)
  cat("Bioconductor:", as.character(BiocManager::version()), "\n")
  for (p in c("DESeq2", "limma", "TCGAbiolinks", "survival", "glmnet", "GEOquery"))
    cat(p, as.character(packageVersion(p)), "\n")
})

safe("TCGA RNA-seq SE (tcga_lihc_rnaseq_se.rds)", {
  se <- readRDS("data/raw/tcga_lihc_rnaseq_se.rds")
  print(dim(se)); print(assayNames(se))
  print(colnames(colData(se)))
  print(head(colnames(se), 3))
  print(table(se$sample_type, useNA = "ifany"))
  cat("Patients with >1 sample:", sum(duplicated(substr(colnames(se), 1, 12))), "\n")
})

safe("Raw RNA-seq per-sample files", {
  f <- list.files("data/raw/TCGA-LIHC/Transcriptome_Profiling", pattern = "star_gene_counts.tsv$",
                  recursive = TRUE, full.names = TRUE)
  cat("n files:", length(f), "\n")
  print(read.delim(f[1], comment.char = "#", nrows = 8, check.names = FALSE))
})

safe("Raw methylation per-sample files", {
  f <- list.files("data/raw/GDCdata/TCGA-LIHC/DNA_Methylation", pattern = "level3betas.txt$",
                  recursive = TRUE, full.names = TRUE)
  cat("n files:", length(f), "\n")
  print(readLines(f[1], 3))
})

safe("Methylation SE (tcga_lihc_meth450_se.rds)", {
  se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
  print(dim(se)); print(assayNames(se))
  print(colnames(colData(se)))
  print(head(colnames(se), 3))
  print(table(se$sample_type, useNA = "ifany"))
})

safe("Methylation filtered (tcga_lihc_meth_filtered.rds)", {
  b <- readRDS("data/raw/tcga_lihc_meth_filtered.rds"); print(class(b)); print(dim(b))
})

safe("Clinical (tcga_lihc_clinical.rds)", {
  cl <- readRDS("data/raw/tcga_lihc_clinical.rds")
  print(class(cl)); print(dim(cl)); print(names(cl))
})

safe("dds_filtered: design, samples, genes", {
  d <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")
  print(dim(d)); print(design(d)); print(colnames(colData(d)))
  print(table(d$sample_type, useNA = "ifany"))
})

safe("dds_DEseqed: result names", {
  d <- readRDS("data/raw/tcga_lihc_dds_DEseqed.rds"); print(resultsNames(d))
})

safe("vsd", {
  v <- readRDS("data/raw/tcga_lihc_vsd.rds"); print(dim(v))
  v2 <- readRDS("data/raw/tcga_lihc_vsd_tumor_matched.rds"); print(dim(v2))
})

safe("TCGA DEG table (both possible locations)", {
  for (p in c("03_results/dge/dge_tcga_tumor_vs_normal.rds", "03_results/dge_tcga_tumor_vs_normal.rds")) {
    cat(p, "exists:", file.exists(p), "\n")
    if (file.exists(p)) { x <- readRDS(p); print(dim(x)); print(names(x)); print(head(x, 3)) }
  }
})

safe("26 nexus gene list", {
  for (p in c("03_results/dge/dge_nexus_genes.csv", "03_results/dge_nexus_genes.csv")) {
    if (file.exists(p)) {
      x <- read.csv(p); cat(p, "rows:", nrow(x), "\n"); print(names(x))
      ch <- names(x)[sapply(x, is.character)][1]; dput(x[[ch]]); break
    }
  }
})

safe("GEO objects", {
  for (g in c("GSE14520", "GSE36376", "GSE76427")) {
    o <- readRDS(sprintf("data/raw/geo/%s.rds", g))
    cat("\n", g, "- ExpressionSets:", length(o), "\n")
    for (i in seq_along(o)) cat("  set", i, "platform", Biobase::annotation(o[[i]]), "dim", dim(o[[i]]), "\n")
    print(head(colnames(Biobase::pData(o[[1]])), 40))
  }
})

safe("Expected paths used by scripts", {
  chk <- c("02_scripts/data/raw/tcga_lihc_maf.rds", "data/raw/tcga_lihc_maf.rds",
           "03_results/dmeth_nexus_summary.csv",
           "03_results/association/17_GSE56588_DUSP10_methylation.csv",
           "03_results/association/19_DUSP10_probe_table.csv",
           "03_results/methylation/probes_rho_beyond_0.3.csv")
  print(data.frame(file = chk, exists = file.exists(chk)))
  cat("\nFiles in 03_results/association:\n"); print(list.files("03_results/association"))
  cat("\nFiles in 03_results/dge:\n"); print(list.files("03_results/dge"))
})

sink()
cat("\nDone. Output saved to 05_validation/S0_probe_output.txt\n")