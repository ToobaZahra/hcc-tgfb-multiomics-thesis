## 64: Boxplots, tumor vs normal, TCGA-LIHC, for DUSP10, SMAD2, MYC, MAPK8 and HDAC11
##   Figure 1: expression (VST)
##   Figure 2: gene-level promoter methylation (mean beta over the gene's probes)
##   Figure 3: methylation at the strongest probe of each gene (DUSP10 = cg23932332)
## Saves PNGs to 04_figures and a table of the p-values and differences to 03_results/association/34_boxplot_stats_5genes.csv
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(ggplot2); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
fig_dir <- "04_figures"; a <- "03_results/association"
dir.create(fig_dir, showWarnings = FALSE)
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")      # edit this list to add genes (e.g. "NPC1")
cols <- c(Normal = "#2C7BB6", Tumor = "#D7301F")

# ---- 1. Expression: tumors (primary, one per patient) vs solid tissue normals ----
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
bcE <- substr(colnames(expr), 14, 15)
et <- which(bcE == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
en <- which(bcE == "11")
ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw <- match(ens, sub("\\..*$", "", rownames(expr)))
if (any(is.na(rw))) stop("Not found in the expression matrix: ", paste(genes[is.na(rw)], collapse = ", "))
ex_long <- bind_rows(lapply(seq_along(genes), function(i) rbind(
  data.frame(gene = genes[i], group = "Normal", value = as.numeric(expr[rw[i], en])),
  data.frame(gene = genes[i], group = "Tumor",  value = as.numeric(expr[rw[i], et])))))

# ---- 2. Methylation: gene-level mean beta and strongest probe ----
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
nor <- which(bc == "11")
pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))
gene_rows <- list(); probe_rows <- list(); key_of <- character(0)
for (G in genes) {
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No methylation probes for", G, ": skipped in the methylation figures.\n"); next }
  key <- if (G == "DUSP10" && "cg23932332" %in% pg$probe) "cg23932332" else pg$probe[which.max(abs(pg$rho))]
  key_of[G] <- key
  gm <- colMeans(assay(se[pg$probe, c(tum, nor)]), na.rm = TRUE)
  km <- as.numeric(assay(se[key, c(tum, nor)]))
  grp <- c(rep("Tumor", length(tum)), rep("Normal", length(nor)))
  gene_rows[[G]] <- data.frame(gene = G, group = grp, value = unname(gm))
  probe_rows[[G]] <- data.frame(gene = paste0(G, "\n", key), group = grp, value = km)
  cat(G, ": ", nrow(pg), " probes; probe shown = ", key, "\n", sep = "")
}
me_gene <- bind_rows(gene_rows); me_probe <- bind_rows(probe_rows)

# ---- 3. Plot helper ----
stats_for <- function(df, measure) {
  df %>% filter(is.finite(value)) %>% group_by(gene) %>% summarise(
    measure = measure,
    n_tumor = sum(group == "Tumor"), n_normal = sum(group == "Normal"),
    mean_tumor = mean(value[group == "Tumor"]), mean_normal = mean(value[group == "Normal"]),
    difference = mean_tumor - mean_normal,
    p_wilcox = wilcox.test(value[group == "Tumor"], value[group == "Normal"])$p.value, .groups = "drop") }

make_plot <- function(df, levs, ylab, title, file, width) {
  df$gene <- factor(df$gene, levs); df$group <- factor(df$group, c("Normal", "Tumor"))
  st <- stats_for(df, "x") %>%
    mutate(label = sprintf("p = %.2g\nn = %d T / %d N", p_wilcox, n_tumor, n_normal), group = "Tumor")
  st$gene <- factor(st$gene, levs)
  p <- ggplot(df[is.finite(df$value), ], aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.3, color = "grey20") +
    geom_text(data = st, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.2, size = 3) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(title = title, subtitle = "TCGA-LIHC; Wilcoxon test (raw p); T = tumor, N = solid tissue normal", x = NULL, y = ylab) +
    theme_minimal(base_size = 11) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, file), p, width = width, height = 4.6, dpi = 300, bg = "white")
  invisible(p) }

make_plot(ex_long, genes, "Expression (VST)", "Expression, tumor vs normal", "fig_5genes_expression_tumor_vs_normal.png", 14)
gl <- intersect(genes, names(gene_rows))
make_plot(me_gene, gl, "Mean beta (promoter probes)", "Gene-level methylation, tumor vs normal",
          "fig_5genes_methylation_gene_level_tumor_vs_normal.png", 14)
make_plot(me_probe, paste0(gl, "\n", key_of[gl]), "Beta value", "Methylation at the strongest probe, tumor vs normal",
          "fig_5genes_methylation_key_probe_tumor_vs_normal.png", 14)

# ---- 4. Table behind the figures ----
tab <- bind_rows(stats_for(ex_long, "Expression (VST)"),
                 stats_for(me_gene, "Gene-level methylation (beta)"),
                 stats_for(me_probe %>% mutate(gene = sub("\n.*", "", gene)), "Strongest-probe methylation (beta)"))
write.csv(tab, file.path(a, "34_boxplot_stats_5genes.csv"), row.names = FALSE)
cat("\n== Values behind the boxplots ==\n"); print(as.data.frame(tab), digits = 3)
cat("\nSaved 3 figures to", fig_dir, "and 34_boxplot_stats_5genes.csv\n")
cat("Raw Wilcoxon p-values; the strongest probe is chosen on TCGA expression correlation, so its panel is descriptive.\n")