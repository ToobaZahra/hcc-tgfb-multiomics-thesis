## 43: separate tables -- expression and methylation vs grade (G1-G4) and vs stage (I-IV)
## gene level (expression, methylation) and probe level (methylation), mean value + n for every grade/stage
library(SummarizedExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
tab_dir <- "03_results/association"; dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── 1. Promoter probes of the 26 genes ───────────────────
se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
ann <- as.data.frame(rowData(se))[, c("gene", "chrm_A", "beg_A")]
ann$probe <- rownames(ann)
g_list <- strsplit(as.character(ann$gene), ";")
hit  <- vapply(g_list, function(g) any(g %in% genes), logical(1))
keep <- hit & !is.na(ann$chrm_A) & !is.na(ann$beg_A)
cand <- ann[keep, ]; cand$gene_list <- g_list[keep]

ent <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")); ent <- ent[!is.na(ent)]
tx_by_gene <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")[ent]
prom <- unlist(promoters(tx_by_gene, upstream = 1500, downstream = 500))
prom$symbol <- names(ent)[match(names(prom), ent)]
prom <- keepStandardChromosomes(prom, pruning.mode = "coarse")
probe_gr <- GRanges(cand$chrm_A, IRanges(cand$beg_A, cand$beg_A + 1))
ov <- suppressWarnings(findOverlaps(probe_gr, prom))
probe_gene <- unique(data.frame(probe = cand$probe[queryHits(ov)], gene = prom$symbol[subjectHits(ov)]))
probe_gene <- probe_gene[mapply(function(p, g) g %in% cand$gene_list[match(p, cand$probe)][[1]],
                                probe_gene$probe, probe_gene$gene), ]

# ── 2. Tumors, grade (1-4) and stage (1-4) ───────────────
bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cd  <- as.data.frame(colData(se))[tum, ]; cd$patient <- substr(colnames(se)[tum], 1, 12)
grade_num <- match(cd$tumor_grade, c("G1", "G2", "G3", "G4"))
st <- cd$ajcc_pathologic_stage
stage_num <- ifelse(grepl("^Stage I[ABC]?$", st), 1, ifelse(grepl("^Stage II[ABC]?$", st), 2,
                                                            ifelse(grepl("^Stage III[ABC]?$", st), 3, ifelse(grepl("^Stage IV[ABC]?$", st), 4, NA))))
names(grade_num) <- names(stage_num) <- cd$patient
cat("Grade:\n"); print(table(grade_num, useNA = "ifany"))
cat("Stage (1=I ... 4=IV):\n"); print(table(stage_num, useNA = "ifany"))

# ── 3. Matrices ──────────────────────────────────────────
# probe beta (tumors)
bt <- assay(se[unique(probe_gene$probe), tum]); colnames(bt) <- cd$patient
# gene-level beta = mean of promoter probes
gb <- t(sapply(unique(probe_gene$gene), function(gn)
  colMeans(bt[probe_gene$probe[probe_gene$gene == gn], , drop = FALSE], na.rm = TRUE)))
# gene expression (VST, tumors)
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
idx <- match(sym2ens, sub("\\..*$", "", rownames(expr))); ok <- !is.na(idx)
ex <- expr[idx[ok], , drop = FALSE]; rownames(ex) <- names(sym2ens)[ok]
et <- which(substr(colnames(ex), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(ex)[et], 1, 12))]
ex <- ex[, et, drop = FALSE]; colnames(ex) <- substr(colnames(ex), 1, 12)

# ── 4. Generic table: mean + n per level + Spearman ──────
by_level <- function(M, v, labs, what) {
  v <- v[colnames(M)]
  out <- do.call(rbind, lapply(rownames(M), function(f) {
    x <- M[f, ]; k <- is.finite(x) & !is.na(v)
    ct <- if (sum(k) >= 10 && sd(x[k]) > 0)
      suppressWarnings(cor.test(x[k], v[k], method = "spearman", exact = FALSE)) else NULL
    d <- data.frame(feature = f, stringsAsFactors = FALSE)
    for (i in seq_along(labs)) d[[paste0(what, "_", labs[i])]] <- mean(x[k & v == i])
    for (i in seq_along(labs)) d[[paste0("n_", labs[i])]] <- sum(k & v == i)
    d$rho <- if (is.null(ct)) NA_real_ else unname(ct$estimate)
    d$p   <- if (is.null(ct)) NA_real_ else ct$p.value
    d
  }))
  out$padj <- p.adjust(out$p, "BH")
  out$significant <- !is.na(out$padj) & out$padj < 0.05
  out
}
G <- c("G1", "G2", "G3", "G4"); S <- c("StageI", "StageII", "StageIII", "StageIV")

t1 <- by_level(ex, grade_num, G, "expr") %>% dplyr::rename(gene = feature) %>% dplyr::arrange(padj)
t2 <- by_level(ex, stage_num, S, "expr") %>% dplyr::rename(gene = feature) %>% dplyr::arrange(padj)
t3 <- by_level(gb, grade_num, G, "beta") %>% dplyr::rename(gene = feature) %>% dplyr::arrange(padj)
t4 <- by_level(gb, stage_num, S, "beta") %>% dplyr::rename(gene = feature) %>% dplyr::arrange(padj)
t5 <- by_level(bt, grade_num, G, "beta") %>% dplyr::rename(probe = feature) %>%
  dplyr::left_join(probe_gene, by = "probe") %>% dplyr::relocate(gene, .before = probe) %>% dplyr::arrange(padj)
t6 <- by_level(bt, stage_num, S, "beta") %>% dplyr::rename(probe = feature) %>%
  dplyr::left_join(probe_gene, by = "probe") %>% dplyr::relocate(gene, .before = probe) %>% dplyr::arrange(padj)

write.csv(t1, file.path(tab_dir, "1_expr_vs_grade_gene.csv"),  row.names = FALSE)
write.csv(t2, file.path(tab_dir, "2_expr_vs_stage_gene.csv"),  row.names = FALSE)
write.csv(t3, file.path(tab_dir, "3_meth_vs_grade_gene.csv"),  row.names = FALSE)
write.csv(t4, file.path(tab_dir, "4_meth_vs_stage_gene.csv"),  row.names = FALSE)
write.csv(t5, file.path(tab_dir, "5_meth_vs_grade_probe.csv"), row.names = FALSE)
write.csv(t6, file.path(tab_dir, "6_meth_vs_stage_probe.csv"), row.names = FALSE)

cat("\n1. Expression vs grade (gene)\n");  print(as.data.frame(t1), digits = 3)
cat("\n2. Expression vs stage (gene)\n");  print(as.data.frame(t2), digits = 3)
cat("\n3. Methylation vs grade (gene)\n"); print(as.data.frame(t3), digits = 3)
cat("\n4. Methylation vs stage (gene)\n"); print(as.data.frame(t4), digits = 3)
cat("\n5. Methylation vs grade (probe):", nrow(t5), "probes,", sum(t5$significant), "significant\n")
cat("6. Methylation vs stage (probe):", nrow(t6), "probes,", sum(t6$significant), "significant\n")
cat("Tables saved in", tab_dir, "\n")