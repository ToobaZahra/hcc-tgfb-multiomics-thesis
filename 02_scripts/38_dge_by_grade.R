## 38: DESeq2 High vs Low grade (all genes) + stage/grade table for all 26 nexus genes
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/dge"

nexus <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Grade per patient (taken from the methylation SE clinical columns) ──
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
cdm <- as.data.frame(colData(se))
cdm$patient <- substr(colnames(se), 1, 12)
cdm <- cdm[substr(colnames(se), 14, 15) == "01", ]        # tumor columns only
cdm <- cdm[!duplicated(cdm$patient), ]
grade_lookup <- setNames(cdm$tumor_grade, cdm$patient)

# ── 2. DESeq2 object ─────────────────────────────────────
dds <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")
dds$patient <- substr(dds$barcode, 1, 12)
dds$grade   <- grade_lookup[dds$patient]
# Low = G1+G2, High = G3+G4 (GX / missing excluded)
dds$grade_group <- ifelse(dds$grade %in% c("G1", "G2"), "Low",
                          ifelse(dds$grade %in% c("G3", "G4"), "High", NA))

keep  <- dds$sample_type == "Primary Tumor" & !is.na(dds$grade_group)
tumor <- dds[, keep]
tumor$grade_group <- factor(tumor$grade_group, levels = c("Low", "High"))
cat("Grade group distribution:\n"); print(table(tumor$grade_group))

# ── 3. DESeq2: High vs Low ───────────────────────────────
design(tumor) <- ~ grade_group
tumor <- DESeq(tumor)
res <- results(tumor, contrast = c("grade_group", "High", "Low"))
res_df <- as.data.frame(res)
res_df$ensembl_id <- gsub("\\..*$", "", rownames(res_df))
res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$ensembl_id,
                        column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
saveRDS(res_df, file.path(out_dir, "dge_high_vs_low_grade.rds"))
write.csv(res_df, file.path(out_dir, "dge_high_vs_low_grade.csv"), row.names = FALSE)
summary(res)

# ── 4. Table for all 26 nexus genes: tumor vs normal, stage, grade ──
best <- function(df, tag) {
  df %>% filter(symbol %in% nexus) %>%
    group_by(symbol) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(gene = symbol, log2FC = log2FoldChange, padj = padj) %>%
    rename_with(~ paste0(.x, "_", tag), c(log2FC, padj))
}
stage <- readRDS(file.path(out_dir, "dge_late_vs_early.rds"))
tab <- best(stage, "stage") %>%
  full_join(best(res_df, "grade"), by = "gene")

tvn_file <- file.path(out_dir, "dge_tcga_tumor_vs_normal.rds")
if (file.exists(tvn_file)) {
  tvn <- readRDS(tvn_file)
  tab <- full_join(best(tvn, "tumor_vs_normal"), tab, by = "gene")
}
tab <- tab %>% arrange(padj_grade)
tab$stage_sig <- !is.na(tab$padj_stage) & tab$padj_stage < 0.05
tab$grade_sig <- !is.na(tab$padj_grade) & tab$padj_grade < 0.05
write.csv(tab, file.path(out_dir, "nexus_stage_grade_dge.csv"), row.names = FALSE)
print(as.data.frame(tab), digits = 3)
cat("\nGenes missing from table:", paste(setdiff(nexus, tab$gene), collapse = ", "), "\n")
