## 52: forest plot of survival HRs (per 1 SD), TCGA-LIHC, GSE14520, GSE76427
library(dplyr); library(ggplot2)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"; fig_dir <- "04_figures"
genes_plot <- c("DUSP10", "NPC1", "SMAD2")

# 1. HR, CI, n, deaths per cohort
tc <- read.csv(file.path(a, "8_expr_vs_survival_gene.csv")) %>%
  transmute(gene, cohort = "TCGA-LIHC", HR, lo = HR_low, hi = HR_high, p = p_cox, n, events)
ge <- read.csv(file.path(a, "15_GEO_test_results_long.csv")) %>%
  filter(!is.na(HR)) %>%
  transmute(gene, cohort = dataset, HR, lo = HR_low, hi = HR_high, p = cox_p, n = n_surv, events)
D <- bind_rows(tc, ge) %>% filter(gene %in% genes_plot)

# 2. Pooled GEO (fixed effect, inverse variance) + Holm + heterogeneity
pool <- D %>% filter(cohort != "TCGA-LIHC") %>% group_by(gene) %>%
  mutate(lh = log(HR), se = (log(hi) - log(lo)) / (2 * 1.96), w = 1 / se^2) %>%
  summarise(lh_p = sum(w * lh) / sum(w), se_p = sqrt(1 / sum(w)),
            Q = sum(w * (lh - sum(w * lh) / sum(w))^2), k = n(),
            n = sum(n), events = sum(events), .groups = "drop") %>%
  mutate(I2 = pmax(0, (Q - (k - 1)) / Q) * 100,
         Q_p = pchisq(Q, df = k - 1, lower.tail = FALSE),
         p = 2 * pnorm(-abs(lh_p / se_p)),
         p_holm = p.adjust(p, method = "holm")) %>%
  transmute(gene, cohort = "GEO pooled", HR = exp(lh_p),
            lo = exp(lh_p - 1.96 * se_p), hi = exp(lh_p + 1.96 * se_p),
            p, n, events, p_holm, I2, Q_p)
D <- bind_rows(D, pool)

lv <- c("TCGA-LIHC", "GSE14520", "GSE76427", "GEO pooled")
D$cohort <- factor(D$cohort, rev(lv))
D$gene   <- factor(D$gene, genes_plot)
D <- D %>% arrange(gene, cohort)

D$label <- sprintf("HR %.2f (%.2f-%.2f), p = %.2g", D$HR, D$lo, D$hi, D$p)
isp <- D$cohort == "GEO pooled"
D$label[isp] <- paste0(D$label[isp],
                       sprintf("\nHolm p = %.2g; I2 = %.0f%% (Q p = %.2g)", D$p_holm[isp], D$I2[isp], D$Q_p[isp]))
D$ylab <- factor(paste0(D$cohort, "\n(n = ", D$n, ", ", D$events, " deaths)"),
                 levels = unique(paste0(D$cohort, "\n(n = ", D$n, ", ", D$events, " deaths)")))
D$type <- ifelse(isp, "pooled", "cohort")

out <- D %>% dplyr::select(gene, cohort, HR, lo, hi, p, n, events, p_holm, I2, Q_p)
write.csv(out, file.path(a, "16_forest_survival_table.csv"), row.names = FALSE)
cat("== Survival hazard ratios (per 1 SD, unadjusted) ==\n")
print(as.data.frame(out), digits = 3)

# 3. Figure
f <- ggplot(D, aes(HR, ylab, color = cohort, shape = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.2) +
  geom_text(aes(label = label, x = 2.3), hjust = 0, size = 3, color = "grey20", lineheight = 0.95) +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2), limits = c(0.4, 7)) +
  scale_color_manual(values = c(`TCGA-LIHC` = "#D7301F", GSE14520 = "#2C7BB6",
                                GSE76427 = "#4DAF4A", `GEO pooled` = "black"), guide = "none") +
  scale_shape_manual(values = c(cohort = 16, pooled = 18), guide = "none") +
  labs(title = "Expression and overall survival across independent cohorts",
       subtitle = paste0("Unadjusted hazard ratio per 1 SD of expression (HR < 1: higher expression, better survival)\n",
                         "GEO pooled = GSE14520 + GSE76427, fixed effect. Cohort p-values nominal; Holm over 3 genes on pooled row"),
       x = "Hazard ratio (log scale)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "fig_forest_survival.png"), f, width = 10, height = 9, dpi = 300, bg = "white")
cat("Saved fig_forest_survival.png and 16_forest_survival_table.csv\n")