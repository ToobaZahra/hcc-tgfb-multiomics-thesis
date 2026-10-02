library(dplyr)
library(tidyr)
library(ggplot2)

# ── Setup ────────────────────────────────────────────────
setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")   # project root
out_dir <- "03_results/variants/somatic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ─────────────────────────────────────────────────
maf <- readRDS("02_scripts/data/raw/tcga_lihc_maf.rds")

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

coding_nonsilent <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                      "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                      "Splice_Site","Translation_Start_Site","Nonstop_Mutation")

# ── Sanity checks ────────────────────────────────────────
maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)
cat("Patients in MAF:", n_distinct(maf$patient), "\n")
print(table(maf$Variant_Classification))

# ── Filter to nexus genes, coding non-silent ─────────────
nexus_maf <- maf %>%
  filter(Hugo_Symbol %in% nexus_genes,
         Variant_Classification %in% coding_nonsilent) %>%
  select(patient, Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
         Variant_Type, HGVSp_Short, IMPACT, SIFT, PolyPhen, t_depth, t_alt_count)

write.csv(nexus_maf,
          file.path(out_dir, "nexus_somatic_mutations_per_patient.csv"),
          row.names = FALSE)

# ── Patient x gene matrix (all MAF patients, zeros included) ─
all_patients <- data.frame(patient = sort(unique(maf$patient)))

mut_matrix <- nexus_maf %>%
  distinct(patient, Hugo_Symbol) %>%
  mutate(mutated = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mutated, values_fill = 0)

mut_matrix_full <- all_patients %>%
  left_join(mut_matrix, by = "patient") %>%
  mutate(across(-patient, ~ replace_na(.x, 0)))

# make sure all 26 genes have a column, even with zero mutations
missing_genes <- setdiff(nexus_genes, colnames(mut_matrix_full))
for (g in missing_genes) mut_matrix_full[[g]] <- 0
mut_matrix_full <- mut_matrix_full[, c("patient", nexus_genes)]

write.csv(mut_matrix_full,
          file.path(out_dir, "nexus_mutation_matrix.csv"),
          row.names = FALSE)

# ── Summary ──────────────────────────────────────────────
cat("\nMutations per gene:\n")
print(sort(table(factor(nexus_maf$Hugo_Symbol, levels = nexus_genes)),
           decreasing = TRUE))
cat("\nPatients with >=1 nexus mutation:",
    sum(rowSums(mut_matrix_full[, nexus_genes]) > 0),
    "of", nrow(mut_matrix_full), "\n")

# ══ Figures ══════════════════════════════════════════════
n_total <- nrow(mut_matrix_full)

# ── Figure 1: mutation frequency per gene ────────────────
freq <- data.frame(gene = nexus_genes,
                   n = colSums(mut_matrix_full[, nexus_genes])) %>%
  mutate(pct = 100 * n / n_total,
         gene = reorder(gene, n))

p1 <- ggplot(freq, aes(gene, n)) +
  geom_col(fill = "#2C7BB6", width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)),
            hjust = -0.1, size = 3.2) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Somatic mutation frequency in TGF-β nexus genes",
       subtitle = paste0("TCGA-LIHC, n = ", n_total,
                         " patients (coding non-silent mutations)"),
       x = NULL, y = "Patients with mutation") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(out_dir, "fig1_nexus_mutation_frequency.png"),
       p1, width = 7, height = 6, dpi = 300, bg = "white")

# ── Figure 2: oncoplot ───────────────────────────────────
cell <- nexus_maf %>%
  group_by(patient, Hugo_Symbol) %>%
  summarise(type = if (n() > 1) "Multi_hit" else first(Variant_Classification),
            .groups = "drop")

gene_order <- as.character(freq$gene[order(-freq$n)])   # most mutated first

bin <- mut_matrix_full %>% filter(patient %in% cell$patient)
ord <- do.call(order, lapply(gene_order, function(g) -bin[[g]]))
patient_levels <- bin$patient[ord]

bg <- expand.grid(patient = patient_levels, gene = gene_order,
                  stringsAsFactors = FALSE)

cell <- cell %>%
  mutate(patient = factor(patient, levels = patient_levels),
         gene = factor(Hugo_Symbol, levels = rev(gene_order)))
bg$patient <- factor(bg$patient, levels = patient_levels)
bg$gene <- factor(bg$gene, levels = rev(gene_order))

cols <- c(Missense_Mutation = "#33A02C", Nonsense_Mutation = "#E31A1C",
          Frame_Shift_Del = "#1F78B4", Frame_Shift_Ins = "#6A3D9A",
          In_Frame_Del = "#FF7F00", In_Frame_Ins = "#FDBF6F",
          Splice_Site = "#B15928", Translation_Start_Site = "#A6CEE3",
          Nonstop_Mutation = "#FB9A99", Multi_hit = "black")

p2 <- ggplot() +
  geom_tile(data = bg, aes(patient, gene), fill = "grey92", color = "white") +
  geom_tile(data = cell, aes(patient, gene, fill = type), color = "white") +
  scale_fill_manual(values = cols, drop = TRUE, name = "Variant type") +
  labs(title = "Oncoplot: nexus gene mutations in TCGA-LIHC",
       subtitle = paste0(length(patient_levels), " of ", n_total,
                         " patients carry ≥1 mutation"),
       x = paste0(length(patient_levels), " patients"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

ggsave(file.path(out_dir, "fig2_nexus_oncoplot.png"),
       p2, width = 11, height = 6, dpi = 300, bg = "white")

print(p1); print(p2)