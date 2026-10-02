library(SummarizedExperiment); library(dplyr)
setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")

se  <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
sel <- read.csv("03_results/methylation/probes_rho_beyond_0.3.csv")

bc  <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]

nor <- which(substr(colnames(se), 14, 15) == "11")
pr  <- unique(sel$probe)
b   <- assay(se[pr, c(nor, tum)]); nn <- length(nor)

out <- do.call(rbind, lapply(pr, function(p) {
  x <- b[p, seq_len(nn)]; y <- b[p, -seq_len(nn)]
  data.frame(probe = p, beta_normal = mean(x, na.rm = TRUE),
             beta_tumor = mean(y, na.rm = TRUE),
             p = suppressWarnings(wilcox.test(y, x)$p.value))
}))
out$delta <- out$beta_tumor - out$beta_normal
out$padj  <- p.adjust(out$p, "BH")
out$tumor_vs_normal <- case_when(out$padj < 0.05 & out$delta >=  0.10 ~ "Hyper",
                                 out$padj < 0.05 & out$delta <= -0.10 ~ "Hypo",
                                 TRUE ~ "No change")
out$level_in_tumor <- case_when(out$beta_tumor >= 0.6 ~ "High",
                                out$beta_tumor <= 0.2 ~ "Low", TRUE ~ "Intermediate")

final <- sel %>% select(gene, probe, dist_TSS, rho) %>%
  left_join(out %>% select(probe, beta_normal, beta_tumor, delta, padj, tumor_vs_normal, level_in_tumor),
            by = "probe") %>% arrange(rho)
write.csv(final, "03_results/methylation/probes_rho_beyond_0.3_methylation.csv", row.names = FALSE)
print(as.data.frame(final))
cat("Normal samples:", nn, "\n")