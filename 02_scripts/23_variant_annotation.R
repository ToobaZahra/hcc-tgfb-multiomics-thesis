# 23_variant_annotation.R
library(httr)
library(jsonlite)
library(dplyr)

# Nexus genes
nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# Ensembl REST API — get gene IDs
server <- "https://rest.ensembl.org"

get_ensembl_id <- function(gene_symbol) {
  url <- paste0(server, "/xrefs/symbol/homo_sapiens/", gene_symbol,
                "?content-type=application/json")
  r <- GET(url)
  if (status_code(r) != 200) return(NA)
  data <- fromJSON(content(r, "text", encoding = "UTF-8"))
  if (length(data) == 0) return(NA)
  return(data$id[1])
}

cat("Fetching Ensembl IDs for nexus genes...\n")
gene_ids <- data.frame(
  symbol = nexus_genes,
  ensembl_id = sapply(nexus_genes, get_ensembl_id)
)

cat("Done.\n")
print(gene_ids)
write.csv(gene_ids, "03_results/nexus_ensembl_ids.csv", row.names = FALSE)


get_variants <- function(ensembl_id, symbol) {
  url <- paste0(server, "/overlap/id/", ensembl_id,
                "?feature=variation;content-type=application/json")
  r <- tryCatch(
    GET(url, timeout(120)),
    error = function(e) return(NULL)
  )
  if (is.null(r) || status_code(r) != 200) return(NULL)
  data <- tryCatch(
    fromJSON(content(r, "text", encoding = "UTF-8")),
    error = function(e) return(NULL)
  )
  if (is.null(data) || length(data) == 0) return(NULL)
  data$gene_symbol <- symbol
  return(data)
}

cat("Fetching variants for nexus genes...\n")
all_variants <- list()

for (i in 1:nrow(gene_ids)) {
  sym <- gene_ids$symbol[i]
  ens <- gene_ids$ensembl_id[i]
  cat("Fetching:", sym, "\n")
  vars <- get_variants(ens, sym)
  if (!is.null(vars)) {
    all_variants[[sym]] <- vars
    cat(" →", nrow(vars), "variants\n")
  } else {
    cat(" → failed or empty\n")
  }
  Sys.sleep(1)
}

variants_df <- bind_rows(all_variants)
cat("Total variants fetched:", nrow(variants_df), "\n")

saveRDS(variants_df, "03_results/nexus_variants_raw.rds")
cat("Variants saved.\n")


# Save RDS (already done successfully)
# For CSV, flatten list columns first
variants_flat <- variants_df
list_cols <- sapply(variants_flat, is.list)
variants_flat[list_cols] <- lapply(variants_flat[list_cols], 
                                   function(x) sapply(x, paste, collapse = ";"))

write.csv(variants_flat, "03_results/nexus_variants_raw.csv", row.names = FALSE)
cat("Variants CSV saved.\n")
cat("Columns:", paste(colnames(variants_df), collapse = ", "), "\n")


# Filter to clinically relevant variants
cat("Consequence types available:\n")
print(table(variants_df$consequence_type))

cat("\nClinical significance available:\n")
print(table(unlist(variants_df$clinical_significance)))




# Prioritize high-impact variants
high_impact <- c("missense_variant", "frameshift_variant", "stop_gained",
                 "splice_acceptor_variant", "splice_donor_variant",
                 "inframe_deletion", "inframe_insertion", "start_lost", "stop_lost")

pathogenic_sig <- c("pathogenic", "likely pathogenic")

# Filter
variants_priority <- variants_df %>%
  filter(consequence_type %in% high_impact) %>%
  mutate(is_pathogenic = sapply(clinical_significance, function(x) 
    any(x %in% pathogenic_sig))) %>%
  arrange(desc(is_pathogenic), gene_symbol)

cat("High-impact variants:", nrow(variants_priority), "\n")
cat("Pathogenic/likely pathogenic:", sum(variants_priority$is_pathogenic), "\n")
cat("\nPer gene:\n")
print(table(variants_priority$gene_symbol))

# Save
variants_priority_flat <- variants_priority
list_cols <- sapply(variants_priority_flat, is.list)
variants_priority_flat[list_cols] <- lapply(variants_priority_flat[list_cols],
                                            function(x) sapply(x, paste, collapse = ";"))

write.csv(variants_priority_flat, "03_results/nexus_variants_prioritized.csv", 
          row.names = FALSE)
cat("Prioritized variants saved.\n")


# ── Pathogenic variant summary per gene ──────────────────
library(dplyr)

pathogenic_summary <- variants_priority %>%
  group_by(gene_symbol) %>%
  summarise(
    total_high_impact = n(),
    n_pathogenic = sum(is_pathogenic),
    n_missense = sum(consequence_type == "missense_variant"),
    n_frameshift = sum(consequence_type == "frameshift_variant"),
    n_stop_gained = sum(consequence_type == "stop_gained"),
    n_splice = sum(consequence_type %in% c("splice_acceptor_variant",
                                           "splice_donor_variant"))
  ) %>%
  arrange(desc(n_pathogenic))

write.csv(pathogenic_summary, "03_results/nexus_variants_summary.csv", 
          row.names = FALSE)

cat("Pathogenic variant summary:\n")
print(pathogenic_summary, n = 26)