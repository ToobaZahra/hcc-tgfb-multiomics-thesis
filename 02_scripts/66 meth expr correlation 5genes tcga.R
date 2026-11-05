## 66: Methylation vs expression correlation (scatter + Spearman) for DUSP10, SMAD2, MYC, MAPK8, HDAC11
##   TCGA-LIHC only (the only cohort with both layers on the same patients), primary tumors only,
##   one sample per patient, expression = VST, methylation = beta.
## Figures: gene-level mean beta, and the strongest-probe version (same probe choice as script 65).
## Saves PNGs to 04_figures and the numbers to 03_results/association/36_meth_expr_corr_stats_tcga.csv
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(ggplot2); library(dplyr) })

# Save to: E:/MS-Project/hcc-tgfb-multiomics-thesis/02_scripts/66_meth_expr_correlation_5genes_tcga.R
setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
fig_dir <- "04_figures"; a <- "03_results/association"
dir.create(fig_dir, showWarnings = FALSE)
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")

# ---- expression (VST), primary tumors, one per patient ----
vsd  <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
bcE  <- substr(colnames(expr), 14, 15)
et   <- which(bcE == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
pat_e <- substr(colnames(expr)[et], 1, 12)
ens  <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw   <- match(ens, sub("\\..*$", "", rownames(expr)))
if (any(is.na(rw))) stop("Not found in the TCGA expression matrix: ", paste(genes[is.na(rw)], collapse = ", "))

# ---- methylation, primary tumors, one per patient ----
se   <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc   <- substr(colnames(se), 14, 15)
tum  <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
pat_m <- substr(colnames(se)[tum], 1, 12)

common <- intersect(pat_e, pat_m)
ie <- et[match(common, pat_e)]; im <- tum[match(common, pat_m)]
cat("Patients with both expression and methylation (tumor):", length(common), "\n")

# ---- probes per gene (from the probe-level correlation table used in script 65) ----
pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))
gprobes <- list(); keyp <- character(0)
for (G in genes) {
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No methylation probes for", G, ": skipped.\n"); next }
  gprobes[[G]] <- pg$probe
  keyp[G] <- if (G == "DUSP10" && "cg23932332" %in% pg$probe) "cg23932332" else pg$probe[which.max(abs(pg$rho))]
}
mg <- names(gprobes)
beta <- assay(se)

build <- function(use_probe) {
  bind_rows(lapply(mg, function(G) {
    probes <- if (use_probe) keyp[[G]] else gprobes[[G]]
    b <- if (length(probes) == 1) as.numeric(beta[probes, im]) else colMeans(beta[probes, im, drop = FALSE], na.rm = TRUE)
    data.frame(gene = G, beta = b, expr = as.numeric(expr[rw[match(G, genes)], ie]),
               label = if (use_probe) paste0(G, "\n", keyp[[G]]) else G)
  })) %>% filter(is.finite(beta), is.finite(expr))
}

stats_of <- function(df) df %>% group_by(gene, label) %>% summarise(
  n = n(), rho = suppressWarnings(cor(beta, expr, method = "spearman")),
  p_spearman = suppressWarnings(cor.test(beta, expr, method = "spearman", exact = FALSE)$p.value), .groups = "drop")

plot_corr <- function(df, st, title, xlab, file) {
  lv <- unique(df$label); df$label <- factor(df$label, lv); st$label <- factor(st$label, lv)
  st$txt <- sprintf("rho = %.2f\np = %.2g\nn = %d", st$rho, st$p_spearman, st$n)
  p <- ggplot(df, aes(beta, expr)) +
    geom_point(size = 0.9, alpha = 0.45, color = "grey25") +
    geom_smooth(method = "lm", se = TRUE, color = "#D7301F", fill = "#D7301F", alpha = 0.15, linewidth = 0.8) +
    geom_text(data = st, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = -0.1, vjust = 1.2, size = 3) +
    facet_wrap(~label, scales = "free", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = title, subtitle = "TCGA-LIHC primary tumors; Spearman (raw p); line = linear fit",
         x = xlab, y = "Expression (VST)") +
    theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, file), p, width = 10, height = 7, dpi = 300, bg = "white") }

dg <- build(FALSE); sg <- stats_of(dg)
dp <- build(TRUE);  sp <- stats_of(dp)
plot_corr(dg, sg, "Methylation vs expression (gene-level), TCGA-LIHC tumors",
          "Mean beta (gene's probes)", "fig_5genes_meth_vs_expr_gene_level_tcga.png")
plot_corr(dp, sp, "Methylation vs expression (strongest probe), TCGA-LIHC tumors",
          "Beta value (one probe per gene)", "fig_5genes_meth_vs_expr_key_probe_tcga.png")

tab <- bind_rows(sg %>% mutate(measure = "Gene-level mean beta"),
                 sp %>% mutate(measure = "Strongest probe", gene = sub("\n.*", "", gene)))
write.csv(tab, file.path(a, "36_meth_expr_corr_stats_tcga.csv"), row.names = FALSE)
cat("\n== Spearman correlations ==\n"); print(as.data.frame(tab), digits = 3)
cat("\nSaved 2 figures to", fig_dir, "and 36_meth_expr_corr_stats_tcga.csv\n")
cat("Raw p-values. The strongest probe was chosen by correlation in TCGA, so its rho is optimistic.\n")