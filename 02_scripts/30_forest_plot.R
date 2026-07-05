# 30_forest_plot.R
library(ggplot2)

cox_res <- read.csv("03_results/survival/cox_univariate.csv")
sig <- cox_res[cox_res$pval < 0.1, ]

p <- ggplot(sig, aes(x = HR, y = reorder(gene, HR),
                     xmin = CI_low, xmax = CI_high,
                     color = HR > 1)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                     labels = c("Protective", "Risk"),
                     name = "Direction") +
  labs(title = "Forest Plot — Cox Univariate Survival",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_bw()

png("04_figures/forest_plot.png", width = 8, height = 6, units = "in", res = 300)
print(p)
dev.off()
cat("Forest plot saved.\n")