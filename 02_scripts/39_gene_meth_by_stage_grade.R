## 39: gene-level promoter methylation (mean beta of promoter probes) vs stage and grade, all 26 genes
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/methylation"

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Promoter probes (same logic as scripts 33/35) ─────
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit  <- vapply(g_list, function(g) any(g %in% genes), logical(1))
keep <- hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)
cand <- ann[keep, ]; cand$gene_list <- g_list[keep]

ent <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL"))
ent <- ent[!is.na(ent)]
tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")[ent]
prom <- unlist(promoters(tx_by_gene, upstream = 1500, downstream = 500))
prom$symbol <- names(ent)[match(names(prom), ent)]
prom <- keepStandardChromosomes(prom, pruning.mode = "coarse")

probe_gr <- GRanges(cand$chrm_A, IRanges(cand$beg_A, cand$beg_A + 1))
ov <- suppressWarnings(findOverlaps(probe_gr, prom))
probe_gene <- unique(data.frame(probe = cand$probe[queryHits(ov)],
                                gene  = prom$symbol[subjectHits(ov)]))
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]
cat("Promoter probes per gene:\n"); print(table(probe_gene$gene))

# ── 2. Tumor samples + clinical groups ───────────────────
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cd  <- as.data.frame(colData(se))[tum, ]
cd$stage <- ifelse(grepl("^Stage I{1,2}[ABC]?$", cd$ajcc_pathologic_stage), "Early",
                   ifelse(grepl("^Stage (III|IV)", cd$ajcc_pathologic_stage), "Late", NA))
cd$grade <- ifelse(cd$tumor_grade %in% c("G1", "G2"), "Low",
                   ifelse(cd$tumor_grade %in% c("G3", "G4"), "High", NA))

# ── 3. Gene-level beta per tumor patient ─────────────────
bt <- assay(se[unique(probe_gene$probe), tum])
gene_beta <- t(sapply(unique(probe_gene$gene), function(g)
  colMeans(bt[probe_gene$probe[probe_gene$gene == g], , drop = FALSE], na.rm = TRUE)))

# ── 4. Compare groups ────────────────────────────────────
cmp <- function(v, grp, lv) {
  a <- v[which(grp == lv[1])]; b <- v[which(grp == lv[2])]
  data.frame(m1 = mean(a, na.rm = TRUE), m2 = mean(b, na.rm = TRUE),
             p = if (sum(is.finite(a)) >= 3 && sum(is.finite(b)) >= 3)
               suppressWarnings(wilcox.test(a, b)$p.value) else NA_real_)
}
res <- do.call(rbind, lapply(rownames(gene_beta), function(g) {
  s <- cmp(gene_beta[g, ], cd$stage, c("Early", "Late"))
  r <- cmp(gene_beta[g, ], cd$grade, c("Low", "High"))
  data.frame(gene = g, n_probes = sum(probe_gene$gene == g),
             beta_early = s$m1, beta_late = s$m2, delta_stage = s$m2 - s$m1, p_stage = s$p,
             beta_gradeLow = r$m1, beta_gradeHigh = r$m2, delta_grade = r$m2 - r$m1, p_grade = r$p)
}))
res$padj_stage <- p.adjust(res$p_stage, "BH")
res$padj_grade <- p.adjust(res$p_grade, "BH")
res <- res %>% arrange(padj_stage)
res$stage_sig <- res$padj_stage < 0.05
res$grade_sig <- res$padj_grade < 0.05
write.csv(res, file.path(out_dir, "gene_meth_by_stage_grade.csv"), row.names = FALSE)
print(as.data.frame(res %>% dplyr::select(gene, n_probes, delta_stage, padj_stage, delta_grade, padj_grade)), digits = 3)
cat("\nGenes missing:", paste(setdiff(genes, res$gene), collapse = ", "), "\n")