## ============================================================
## Script: 26_gwas_nexus_overlap.R
## Purpose: Overlap HCC GWAS summary stats (GCST90860789/790/791)
##          with the 26 TGF-beta nexus genes, add as a 5th
##          multi-omics layer (GWAS_hit) alongside expression,
##          methylation, survival, and replication.
## Data: GWAS Catalog, Chinaka et al. 2026, Hum Genet Genomics Adv
##       Multi-ancestry meta-analysis of HCC (METAL), GRCh38
## ============================================================

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")

library(data.table)

## --- ONE-TIME SETUP (run once per machine, then comment out) ---
## Offline annotation packages — no network call at runtime, avoids Ensembl outages.
## If any of these are missing you'll get "no package called X" — install in this
## order (rtracklayer is a hidden dependency of GenomicFeatures, install it first
## if you hit that specific error):
##   options(timeout = 600)  # these are large downloads, raise the default 60s timeout
##   BiocManager::install(c("rtracklayer", "GenomicFeatures",
##                           "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
##                        update = FALSE, ask = FALSE)
## -----------------------------------------------------------------

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)

## ------------------------------------------------------------
## Step 1: Get GRCh38 coordinates for the 26 nexus genes
## ------------------------------------------------------------

nexus_genes <- c(
  "TGFB1", "TGFBR1", "TGFBR2",
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  "MAPK8", "MAPK9", "MAPK10",
  "DUSP1", "DUSP4", "DUSP10",
  "MYC", "CDKN1A", "CDKN2B", "SNAI1",
  "TNF", "IL6", "IL10", "IL37",
  "HDAC11", "NPC1", "CCDC110", "TGFBRAP1", "KLF4"
)

# sanity check — should be 26
stopifnot(length(nexus_genes) == 26)

# IL37 is stored as IL1F7 in older annotation builds -- map both, keep whichever resolves
symbol_lookup <- c(nexus_genes, "IL1F7")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_ranges <- genes(txdb)  # GRanges keyed by Entrez ID

# map symbols -> Entrez IDs
entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = symbol_lookup,
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]

matched <- gene_ranges[names(gene_ranges) %in% entrez_map$ENTREZID]
entrez_map_matched <- entrez_map[entrez_map$ENTREZID %in% names(matched), ]

gene_coords <- data.frame(
  hgnc_symbol = entrez_map_matched$SYMBOL[match(names(matched), entrez_map_matched$ENTREZID)],
  chromosome_name = sub("^chr", "", as.character(seqnames(matched))),
  start_position = start(matched),
  end_position = end(matched),
  stringsAsFactors = FALSE
)

# fold IL1F7 rows back to IL37 so downstream code sees the expected symbol
gene_coords$hgnc_symbol[gene_coords$hgnc_symbol == "IL1F7"] <- "IL37"

# keep only standard chromosomes (1-22, X) — drops scaffold/patch duplicates
gene_coords <- gene_coords[gene_coords$chromosome_name %in% c(1:22, "X"), ]
# some symbols may resolve to >1 Entrez ID (readthrough/pseudogene entries) -- keep widest span
gene_coords <- as.data.table(gene_coords)[
  , .(start_position = min(start_position), end_position = max(end_position)),
  by = .(hgnc_symbol, chromosome_name)
]
gene_coords <- as.data.frame(gene_coords)

missing_genes <- setdiff(nexus_genes, gene_coords$hgnc_symbol)
if (length(missing_genes) > 0) {
  warning("No coordinates found for: ", paste(missing_genes, collapse = ", "),
          " -- check gene symbol aliases manually (e.g. via genecards.org) and add coordinates by hand")
}

# TNF sits in the gene-dense MHC region on chr6 and gets dropped by the
# "genes on both strands / multiple sequences" filter in genes(txdb).
# Add its GRCh38 coordinates manually (Ensembl ENSG00000232810).
if ("TNF" %in% missing_genes) {
  tnf_row <- data.frame(
    hgnc_symbol = "TNF",
    chromosome_name = "6",
    start_position = 31575565,
    end_position = 31578336
  )
  gene_coords <- rbind(gene_coords, tnf_row)
}

# should be all 26 now
stopifnot(nrow(gene_coords) == 26)

# pad +/- 100kb to catch cis-regulatory / eQTL variants
FLANK <- 100000
gene_coords$win_start <- pmax(0, gene_coords$start_position - FLANK)
gene_coords$win_end   <- gene_coords$end_position + FLANK

fwrite(gene_coords, "03_results/pathways/nexus_gene_coords_grch38.csv")

## ------------------------------------------------------------
## Step 2: Load GWAS summary stats (multi-ancestry, primary file)
## ------------------------------------------------------------

# GWAS-SSF format columns: chromosome, base_pair_location, effect_allele,
# other_allele, beta, standard_error, effect_allele_frequency, p_value, rsid, n

gwas_multi <- fread("data/raw/GCST90860789.tsv.gz")

# normalize chromosome column to character, strip any "chr" prefix if present
gwas_multi[, chromosome := as.character(chromosome)]
gwas_multi[, chromosome := sub("^chr", "", chromosome)]

## ------------------------------------------------------------
## Step 3: Overlap function -- pull SNPs inside each gene window
## ------------------------------------------------------------

get_gene_hits <- function(gwas_dt, coords_dt, p_threshold = 1e-5) {
  hits_list <- list()
  for (i in seq_len(nrow(coords_dt))) {
    g <- coords_dt[i, ]
    sub <- gwas_dt[chromosome == as.character(g$chromosome_name) &
                     base_pair_location >= g$win_start &
                     base_pair_location <= g$win_end]
    if (nrow(sub) > 0) {
      sub[, gene := g$hgnc_symbol]
      hits_list[[g$hgnc_symbol]] <- sub
    }
  }
  out <- rbindlist(hits_list, fill = TRUE)
  out[p_value < p_threshold]
}

# suggestive (p < 1e-5) and genome-wide significant (p < 5e-8) tables
nexus_hits_suggestive <- get_gene_hits(gwas_multi, gene_coords, p_threshold = 1e-5)
nexus_hits_gws        <- nexus_hits_suggestive[p_value < 5e-8]

setorder(nexus_hits_suggestive, gene, p_value)
setorder(nexus_hits_gws, gene, p_value)

fwrite(nexus_hits_suggestive, "03_results/pathways/nexus_gwas_hits_suggestive_multiancestry.csv")
fwrite(nexus_hits_gws,        "03_results/pathways/nexus_gwas_hits_gws_multiancestry.csv")

cat("Genes with suggestive GWAS signal (p<1e-5):",
    length(unique(nexus_hits_suggestive$gene)), "/ 26\n")
cat("Genes with genome-wide significant signal (p<5e-8):",
    length(unique(nexus_hits_gws$gene)), "/ 26\n")

## ------------------------------------------------------------
## Step 4: Ancestry cross-check (EUR = 790, EAS = 791)
## ------------------------------------------------------------

gwas_eur <- fread("data/raw/GCST90860790.tsv.gz")
gwas_eas <- fread("data/raw/GCST90860791.tsv.gz")

gwas_eur[, chromosome := sub("^chr", "", as.character(chromosome))]
gwas_eas[, chromosome := sub("^chr", "", as.character(chromosome))]

eur_hits <- get_gene_hits(gwas_eur, gene_coords, p_threshold = 1e-5)
eas_hits <- get_gene_hits(gwas_eas, gene_coords, p_threshold = 1e-5)

# match multi-ancestry top SNP per gene against EUR/EAS by rsid, compare direction
top_snp_per_gene <- nexus_hits_suggestive[, .SD[which.min(p_value)], by = gene]

check_direction <- function(top_dt, other_dt, label) {
  merged <- merge(top_dt[, .(gene, rsid, beta_multi = beta)],
                  other_dt[, .(rsid, beta_other = beta, p_value)],
                  by = "rsid", all.x = TRUE)
  merged[, same_direction := sign(beta_multi) == sign(beta_other)]
  setnames(merged, c("beta_other", "p_value", "same_direction"),
           paste0(c("beta_", "p_", "same_dir_"), label))
  merged
}

replication_eur <- check_direction(top_snp_per_gene, eur_hits, "EUR")
replication_eas <- check_direction(top_snp_per_gene, eas_hits, "EAS")

ancestry_check <- Reduce(function(x, y) merge(x, y, by = c("gene", "rsid", "beta_multi"), all.x = TRUE),
                         list(top_snp_per_gene[, .(gene, rsid, beta_multi = beta)],
                              replication_eur, replication_eas))

fwrite(ancestry_check, "03_results/pathways/nexus_gwas_ancestry_check.csv")

## ------------------------------------------------------------
## Step 5: Build GWAS_hit column for the multi-omics integration table
## ------------------------------------------------------------

gwas_layer <- data.table(gene = nexus_genes)
gwas_layer <- merge(gwas_layer, top_snp_per_gene[, .(gene, rsid, p_value, beta)],
                    by = "gene", all.x = TRUE)
gwas_layer[, GWAS_hit := fifelse(!is.na(p_value) & p_value < 5e-8, "GWS",
                                 fifelse(!is.na(p_value) & p_value < 1e-5, "Suggestive", "None"))]
setnames(gwas_layer, c("rsid", "p_value", "beta"), c("top_snp", "gwas_p", "gwas_beta"))

fwrite(gwas_layer, "03_results/pathways/nexus_gwas_layer_for_integration.csv")

cat("\nDone. Merge 03_results/pathways/nexus_gwas_layer_for_integration.csv\n",
    "into your Step 14 multi-omics scoring table as the 5th layer.\n")

## ------------------------------------------------------------
## Step 6: Bar graph
## ------------------------------------------------------------

library(ggplot2)

gene_minp$gene <- factor(gene_minp$gene, levels = gene_minp$gene[order(gene_minp$neglog10p)])

p <- ggplot(gene_minp, aes(x = gene, y = neglog10p)) +
  geom_col(fill = "#4472C4", width = 0.7) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(1e-5), color = "orange", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 2, y = -log10(5e-8) + 0.3, label = "Genome-wide significant (p=5e-8)",
           color = "red", hjust = 0, size = 3.2) +
  annotate("text", x = 2, y = -log10(1e-5) + 0.3, label = "Suggestive (p=1e-5)",
           color = "orange", hjust = 0, size = 3.2) +
  coord_flip() +
  labs(title = "GWAS Signal at 26 TGF-\u03b2 Nexus Genes (HCC, Multi-Ancestry)",
       subtitle = "Best SNP per gene \u00b1100kb window | 0/26 genome-wide significant",
       x = NULL, y = expression(-log[10](p))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave("04_figures/gwas_nexus_genes_barplot.png", p, width = 9, height = 7, dpi = 300)
print(p)