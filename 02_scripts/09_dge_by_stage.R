# 09_dge_by_stage.R
library(DESeq2)

dds <- readRDS("data/raw/tcga_lihc_dds_filtered.rds")
clin <- readRDS("data/raw/tcga_lihc_clinical.rds")

# Map clinical stage to dds samples
dds$patient <- substr(dds$barcode, 1, 12)

# Build stage lookup
stage_lookup <- setNames(clin$ajcc_pathologic_stage, clin$submitter_id)
dds$stage <- stage_lookup[dds$patient]

# Group: I+II = Early, III+IV = Late
early_stages <- c("Stage I", "Stage II")
late_stages  <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                  "Stage IV", "Stage IVA", "Stage IVB")

dds$stage_group <- ifelse(dds$stage %in% early_stages, "Early",
                          ifelse(dds$stage %in% late_stages,  "Late", NA))

# Subset: tumor only with stage info
keep <- dds$sample_type == "Primary Tumor" & !is.na(dds$stage_group)
tumor <- dds[, keep]
tumor$stage_group <- factor(tumor$stage_group, levels = c("Early", "Late"))

cat("Stage group distribution:\n")
print(table(tumor$stage_group))

# Re-run DESeq with new design
design(tumor) <- ~ stage_group
tumor <- DESeq(tumor)

res_stage <- results(tumor, contrast = c("stage_group", "Late", "Early"))
res_stage_df <- as.data.frame(res_stage)
res_stage_df$ensembl_id <- gsub("\\..*$", "", rownames(res_stage_df))

library(org.Hs.eg.db)
res_stage_df$symbol <- mapIds(org.Hs.eg.db,
                              keys = res_stage_df$ensembl_id,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

saveRDS(res_stage_df, "03_results/dge_late_vs_early.rds")
write.csv(res_stage_df, "03_results/dge_late_vs_early.csv", row.names = FALSE)
summary(res_stage)
cat("Stage-stratified DGE done.\n")