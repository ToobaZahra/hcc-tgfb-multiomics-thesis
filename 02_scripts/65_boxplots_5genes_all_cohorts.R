## 65: Tumor vs normal boxplots for DUSP10, SMAD2, MYC, MAPK8, HDAC11 in ALL cohorts
##   Expression:  TCGA-LIHC (VST), GSE14520 (log2), GSE76427 (log2)
##   Methylation: TCGA-LIHC, GSE56588 (cirrhotic samples excluded), GSE54503 (paired tumor vs adjacent, 60 pairs)
## Figures (rows = genes, columns = cohorts): expression, gene-level methylation, strongest-probe methylation.
## Saves PNGs to 04_figures and the numbers behind them to 03_results/association/35_boxplot_stats_all_cohorts.csv
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(GEOquery); library(ggplot2); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
fig_dir <- "04_figures"; a <- "03_results/association"; g <- "data/raw/geo/"
dir.create(fig_dir, showWarnings = FALSE)
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")        # edit to add genes
cols <- c(Normal = "#2C7BB6", Tumor = "#D7301F")

# one plotting block: values for tumor and normal, with a Wilcoxon p (paired if requested)
piece <- function(cohort, gene, tv, nv, paired = FALSE) {
  if (paired) { ok <- is.finite(tv) & is.finite(nv); tv <- tv[ok]; nv <- nv[ok]
  } else { tv <- tv[is.finite(tv)]; nv <- nv[is.finite(nv)] }
  if (length(tv) < 3 || length(nv) < 3) return(NULL)
  list(df = data.frame(cohort = cohort, gene = gene, group = c(rep("Tumor", length(tv)), rep("Normal", length(nv))),
                       value = c(tv, nv)),
       st = data.frame(cohort = cohort, gene = gene, n_tumor = length(tv), n_normal = length(nv),
                       mean_tumor = mean(tv), mean_normal = mean(nv), difference = mean(tv) - mean(nv),
                       p_wilcox = wilcox.test(tv, nv, paired = paired, exact = FALSE)$p.value, paired = paired)) }
collect <- function(lst) { lst <- Filter(Negate(is.null), lst)
list(df = bind_rows(lapply(lst, `[[`, "df")), st = bind_rows(lapply(lst, `[[`, "st"))) }

# ================= EXPRESSION =================
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
bcE <- substr(colnames(expr), 14, 15)
et <- which(bcE == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
en <- which(bcE == "11")
ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw <- match(ens, sub("\\..*$", "", rownames(expr)))
if (any(is.na(rw))) stop("Not found in the TCGA expression matrix: ", paste(genes[is.na(rw)], collapse = ", "))

# GEO helpers (same as scripts 51 and 62)
prep_annot <- function(f) { tab <- Table(getGEO(filename = f))
sc <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
list(ids = as.character(tab$ID), sy = lapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), trimws)) }
gene_vec <- function(ann, mat, gene) {
  idx <- which(vapply(ann$sy, function(x) gene %in% x, logical(1)))
  hit <- intersect(ann$ids[idx], rownames(mat)); if (length(hit) == 0) return(NULL)
  best <- hit[which.max(rowMeans(mat[hit, , drop = FALSE]))]
  list(x = mat[best, ], probe = best) }
m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
tis <- tolower(trimws(cl$`Tissue Type`))
d1 <- data.frame(id = as.character(cl$Affy_GSM), tumor = ifelse(grepl("non", tis), FALSE, grepl("tumor", tis)))
d1 <- d1[d1$id %in% colnames(m1), ]
ann1 <- prep_annot(paste0(g, "GPL3921.annot.gz")); lin1 <- max(m1, na.rm = TRUE) > 100
m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]]); p2 <- p2[match(colnames(m2), rownames(p2)), ]
tis2 <- tolower(as.character(p2$`tissue:ch1`))
d2 <- data.frame(id = colnames(m2), tumor = ifelse(grepl("non", tis2), FALSE, grepl("tumor", tis2)))
ann2 <- prep_annot(paste0(g, "GPL10558.annot.gz")); lin2 <- max(m2, na.rm = TRUE) > 100
geo_piece <- function(d, m, ann, lin, G, cohort) {
  gv <- gene_vec(ann, m, G); if (is.null(gv)) { cat(G, "not found on the", cohort, "array\n"); return(NULL) }
  x <- unname(gv$x[d$id]); if (lin) x <- log2(pmax(x, 1))
  k <- !is.na(d$tumor)
  piece(cohort, G, x[k & d$tumor], x[k & !d$tumor]) }

ex <- collect(unlist(lapply(seq_along(genes), function(i) list(
  piece("TCGA-LIHC", genes[i], as.numeric(expr[rw[i], et]), as.numeric(expr[rw[i], en])),
  geo_piece(d1, m1, ann1, lin1, genes[i], "GSE14520"),
  geo_piece(d2, m2, ann2, lin2, genes[i], "GSE76427"))), recursive = FALSE))

# ================= METHYLATION =================
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))
gprobes <- list(); keyp <- character(0)
for (G in genes) {
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No methylation probes for", G, ": skipped in the methylation figures.\n"); next }
  gprobes[[G]] <- pg$probe
  keyp[G] <- if (G == "DUSP10" && "cg23932332" %in% pg$probe) "cg23932332" else pg$probe[which.max(abs(pg$rho))]
}
mg <- names(gprobes); need <- unique(c(unlist(gprobes), keyp))
cat("Methylation genes:", paste(mg, collapse = ", "), "| probe shown:", paste(keyp, collapse = ", "), "\n")

# GSE56588 beta values (large file; only the needed probes are kept) and sample groups
B56 <- NULL; grp56 <- NULL
tryCatch({
  dest <- paste0(g, "GSE56588/GSE56588_Heptromics_avrg_beta.txt.gz")
  e56 <- getGEO("GSE56588", GSEMatrix = TRUE, getGPL = FALSE, destdir = g)[[1]]; p56 <- pData(e56)
  grp56 <- data.frame(sample = p56$title, group = dplyr::case_when(grepl("HCC", p56$characteristics_ch1) ~ "Tumor",
                                                                   grepl("normal", p56$characteristics_ch1) ~ "Normal", grepl("Cirrhotic", p56$characteristics_ch1) ~ "Cirrhotic"))
  hdr <- names(data.table::fread(dest, nrows = 0)); keep <- which(!grepl("_pval$", hdr))
  Bf <- data.table::fread(dest, select = keep, data.table = FALSE)
  rownames(Bf) <- Bf[[1]]; B56 <- as.matrix(Bf[, -1]); B56 <- B56[rownames(B56) %in% need, , drop = FALSE]; rm(Bf); gc()
  cat("GSE56588:", ncol(B56), "samples;", table(grp56$group)[c("Tumor", "Normal")], "(tumor, normal)\n")
}, error = function(e) cat("GSE56588 skipped:", conditionMessage(e), "\n"))

# GSE54503 (cached series matrix from script 55): same order-based pairing, verified with genotype probes
B54 <- NULL
f54 <- paste0(g, "GSE54503_series_matrix.txt.gz")
if (file.exists(f54)) {
  e54 <- getGEO(filename = f54, getGPL = FALSE); p54 <- pData(e54); Bx <- exprs(e54)
  chr <- grep(":ch1$", names(p54), value = TRUE)
  txt <- tolower(apply(p54[, c("title", chr), drop = FALSE], 1, paste, collapse = " "))
  p54$tissue <- ifelse(grepl("non-?tumou?r|adjacent|normal|benign|surrounding", txt), "Normal",
                       ifelse(grepl("tumou?r|hcc|carcinoma", txt), "Tumor", NA))
  p54 <- p54[order(rownames(p54)), ]
  ni <- which(p54$tissue == "Normal"); ok <- ni < nrow(p54) & p54$tissue[pmin(ni + 1, nrow(p54))] == "Tumor"
  nor_g <- rownames(p54)[ni[ok]]; tum_g <- rownames(p54)[ni[ok] + 1]
  rs <- grep("^rs", rownames(Bx), value = TRUE)
  same <- sapply(seq_along(nor_g), function(i) cor(Bx[rs, nor_g[i]], Bx[rs, tum_g[i]], use = "complete.obs"))
  other <- sapply(seq_along(nor_g), function(i) cor(Bx[rs, nor_g[i]], Bx[rs, tum_g[(i %% length(tum_g)) + 1]], use = "complete.obs"))
  if (median(same) > median(other) + 0.2) { B54 <- Bx; cat("GSE54503:", length(tum_g), "pairs confirmed by genotype probes.\n")
  } else cat("GSE54503 skipped: genotype check did not confirm the pairing.\n")
} else cat("GSE54503 series matrix not found: GSE54503 skipped.\n")

meth_piece <- function(G, probes, label, which_what) {
  gm <- function(M, cols) if (length(probes) == 1) as.numeric(M[probes, cols]) else colMeans(M[probes, cols, drop = FALSE], na.rm = TRUE)
  out <- list(piece("TCGA-LIHC", label, as.numeric(gm(assay(se), tum)), as.numeric(gm(assay(se), nor))))
  if (!is.null(B56)) { pr <- intersect(probes, rownames(B56))
  if (length(pr) > 0) { tcol <- intersect(grp56$sample[grp56$group == "Tumor"], colnames(B56))
  ncol_ <- intersect(grp56$sample[grp56$group == "Normal"], colnames(B56))
  out <- c(out, list(piece("GSE56588", label, as.numeric(gm(B56[pr, , drop = FALSE], tcol)), as.numeric(gm(B56[pr, , drop = FALSE], ncol_))))) } }
  if (!is.null(B54)) { pr <- intersect(probes, rownames(B54))
  if (length(pr) > 0) out <- c(out, list(piece("GSE54503 (paired)", label,
                                               as.numeric(gm(B54[pr, , drop = FALSE], tum_g)), as.numeric(gm(B54[pr, , drop = FALSE], nor_g)), paired = TRUE))) }
  out }
me_gene  <- collect(unlist(lapply(mg, function(G) meth_piece(G, gprobes[[G]], G)), recursive = FALSE))
me_probe <- collect(unlist(lapply(mg, function(G) meth_piece(G, keyp[[G]], paste0(G, "\n", keyp[[G]]))), recursive = FALSE))

# ================= PLOTS =================
make_grid <- function(obj, row_levels, col_levels, ylab, title, file, height) {
  df <- obj$df; st <- obj$st
  df$gene <- factor(df$gene, row_levels); df$cohort <- factor(df$cohort, col_levels); df$group <- factor(df$group, c("Normal", "Tumor"))
  st$gene <- factor(st$gene, row_levels); st$cohort <- factor(st$cohort, col_levels)
  st$label <- sprintf("p = %.2g\nn = %d T / %d N%s", st$p_wilcox, st$n_tumor, st$n_normal, ifelse(st$paired, " (paired)", ""))
  p <- ggplot(df, aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.3, color = "grey20") +
    geom_text(data = st, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.15, size = 2.8) +
    facet_grid(gene ~ cohort, scales = "free_y") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(title = title, subtitle = "Wilcoxon test (raw p; paired for GSE54503); T = tumor, N = normal or adjacent tissue", x = NULL, y = ylab) +
    theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, file), p, width = 10, height = height, dpi = 300, bg = "white") }

make_grid(ex, genes, c("TCGA-LIHC", "GSE14520", "GSE76427"), "Expression (TCGA: VST; GEO: log2)",
          "Expression, tumor vs normal, all cohorts", "fig_5genes_expression_all_cohorts.png", 12)
mcols <- c("TCGA-LIHC", "GSE56588", "GSE54503 (paired)")
make_grid(me_gene, mg, mcols, "Mean beta (gene's probes)", "Gene-level methylation, tumor vs normal, all cohorts",
          "fig_5genes_methylation_gene_level_all_cohorts.png", 12)
make_grid(me_probe, paste0(mg, "\n", keyp[mg]), mcols, "Beta value (one probe per gene)",
          "Methylation at the strongest probe, tumor vs normal, all cohorts", "fig_5genes_methylation_key_probe_all_cohorts.png", 12)

tab <- bind_rows(ex$st %>% mutate(measure = "Expression"), me_gene$st %>% mutate(measure = "Gene-level methylation"),
                 me_probe$st %>% mutate(measure = "Strongest-probe methylation", gene = sub("\n.*", "", gene)))
write.csv(tab, file.path(a, "35_boxplot_stats_all_cohorts.csv"), row.names = FALSE)
cat("\n== Values behind the boxplots ==\n"); print(as.data.frame(tab), digits = 3)
cat("\nSaved 3 figures to", fig_dir, "and 35_boxplot_stats_all_cohorts.csv\n")
cat("Raw p-values. TCGA expression is VST and GEO is log2: compare direction, not size. The strongest probe was chosen on TCGA.\n")
cat("GSE54503 uses the 60 verified pairs; GSE56588 has about 10 normals; cirrhotic samples are excluded.\n")