# 24_string_network.R
library(httr)
library(jsonlite)
library(dplyr)

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# STRING API
string_api <- "https://string-db.org/api"
species <- 9606  # human

# Get STRING IDs
url <- paste0(string_api, "/json/get_string_ids?identifiers=",
              paste(nexus_genes, collapse = "%0d"),
              "&species=", species)

r <- GET(url)
string_ids <- fromJSON(content(r, "text", encoding = "UTF-8"))
cat("STRING IDs fetched:", nrow(string_ids), "\n")
print(string_ids[, c("queryItem", "stringId", "preferredName")])



# Get interaction network
string_ids_list <- paste(string_ids$stringId, collapse = "%0d")

url_network <- paste0(string_api, "/json/network?identifiers=",
                      string_ids_list,
                      "&species=", species,
                      "&required_score=400")

r_net <- GET(url_network)
network <- fromJSON(content(r_net, "text", encoding = "UTF-8"))
cat("Interactions fetched:", nrow(network), "\n")
cat("Columns:", colnames(network), "\n")

# Save
write.csv(network, "03_results/string_network.csv", row.names = FALSE)
cat("Network saved.\n")






library(igraph)
library(ggraph)
library(ggplot2)

# Build graph
g <- graph_from_data_frame(
  d = network[, c("preferredName_A", "preferredName_B", "score")],
  directed = FALSE
)

# Load multi-omics results for node coloring
integrated <- read.csv("03_results/multiomics_integrated.csv")
node_colors <- data.frame(
  gene = V(g)$name,
  score = integrated$omics_score[match(V(g)$name, integrated$symbol)]
)
node_colors$score[is.na(node_colors$score)] <- 0
V(g)$omics_score <- node_colors$score

# Plot
png("04_figures/string_network.png", width = 12, height = 10, 
    units = "in", res = 300)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = score), alpha = 0.4, color = "grey60") +
  geom_node_point(aes(color = factor(omics_score)), size = 8) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
  scale_color_manual(values = c("0" = "grey80", "1" = "#FEE08B",
                                "2" = "#FC8D59", "3" = "#D73027",
                                "4" = "#7B0000"),
                     name = "Omics score") +
  scale_edge_width(range = c(0.3, 2)) +
  labs(title = "TGF-β Nexus Gene STRING Network",
       subtitle = "Node color = multi-omics score") +
  theme_graph()
dev.off()
cat("Network plot saved.\n")