# 27_hub_genes.R
library(igraph)
library(dplyr)

network <- read.csv("03_results/string_network.csv")

g <- graph_from_data_frame(
  d = network[, c("preferredName_A", "preferredName_B", "score")],
  directed = FALSE
)

# Centrality measures
degree_cent <- degree(g)
betweenness_cent <- betweenness(g, normalized = TRUE)

hub_df <- data.frame(
  gene = names(degree_cent),
  degree = degree_cent,
  betweenness = round(betweenness_cent, 4)
) %>% arrange(desc(degree))

cat("Top 10 hub genes by degree:\n")
print(hub_df[1:10, ])

write.csv(hub_df, "03_results/string_hub_genes.csv", row.names = FALSE)
cat("Hub genes saved.\n")