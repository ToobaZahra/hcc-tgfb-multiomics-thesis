library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/variants/somatic"

maf <- readRDS("02_scripts/data/raw/tcga_lihc_maf.rds")

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

coding_nonsilent <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                      "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                      "Splice_Site","Translation_Start_Site","Nonstop_Mutation")

ann <- maf %>%
  filter(Hugo_Symbol %in% nexus_genes,
         Variant_Classification %in% coding_nonsilent) %>%
  mutate(
    patient   = substr(Tumor_Sample_Barcode, 1, 12),
    VAF       = round(t_alt_count / t_depth, 3),
    dbSNP_RS  = as.character(dbSNP_RS),
    CLIN_SIG  = as.character(CLIN_SIG),
    COSMIC    = as.character(COSMIC),
    gnomAD_AF = suppressWarnings(as.numeric(gnomAD_AF)),
    known_dbSNP   = !is.na(dbSNP_RS) & !dbSNP_RS %in% c("", "novel", "."),
    known_ClinVar = !is.na(CLIN_SIG) & CLIN_SIG != "",
    known_COSMIC  = !is.na(COSMIC) & COSMIC != "",
    pathogenic_ClinVar = grepl("pathogenic", CLIN_SIG, ignore.case = TRUE) &
      !grepl("conflicting|uncertain", CLIN_SIG, ignore.case = TRUE),
    likely_functional = IMPACT == "HIGH" |
      (grepl("deleterious", SIFT) & grepl("damaging", PolyPhen)),
    possible_germline = !is.na(gnomAD_AF) & gnomAD_AF > 0.001
  ) %>%
  select(patient, Hugo_Symbol, Variant_Classification, HGVSp_Short, IMPACT,
         SIFT, PolyPhen, VAF, dbSNP_RS, CLIN_SIG, COSMIC, gnomAD_AF,
         known_dbSNP, known_ClinVar, known_COSMIC, pathogenic_ClinVar,
         likely_functional, possible_germline) %>%
  arrange(desc(likely_functional), Hugo_Symbol)

write.csv(ann, file.path(out_dir, "nexus_somatic_annotated.csv"), row.names = FALSE)

cat("Total somatic nexus mutations:", nrow(ann), "\n")
cat("In dbSNP (has rsID):          ", sum(ann$known_dbSNP), "\n")
cat("Has ClinVar label:            ", sum(ann$known_ClinVar), "\n")
cat("ClinVar pathogenic:           ", sum(ann$pathogenic_ClinVar), "\n")
cat("In COSMIC:                    ", sum(ann$known_COSMIC), "\n")
cat("Likely functional:            ", sum(ann$likely_functional), "\n")
cat("Possible germline (gnomAD>0.1%):", sum(ann$possible_germline), "\n")

print(as.data.frame(ann %>% filter(likely_functional) %>%
                      select(patient, Hugo_Symbol, HGVSp_Short, IMPACT, SIFT, PolyPhen, VAF)))

# figure
library(ggplot2)

# ── Fig 3: functional vs non-functional mutations per gene ──
f3 <- ann %>%
  mutate(class = ifelse(likely_functional, "Predicted functional", "Other coding")) %>%
  count(Hugo_Symbol, class)

gene_ord <- f3 %>% group_by(Hugo_Symbol) %>% summarise(t = sum(n)) %>%
  arrange(t) %>% pull(Hugo_Symbol)
f3$Hugo_Symbol <- factor(f3$Hugo_Symbol, levels = gene_ord)

p3 <- ggplot(f3, aes(Hugo_Symbol, n, fill = class)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Predicted functional" = "#D7301F",
                               "Other coding" = "grey70"), name = NULL) +
  labs(title = "Nexus gene mutations: predicted functional vs other",
       subtitle = "TCGA-LIHC, 50 coding non-silent mutations (mutation events, not patients)",
       x = NULL, y = "Mutations") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

ggsave(file.path(out_dir, "fig3_functional_vs_other.png"),
       p3, width = 7, height = 5.5, dpi = 300, bg = "white")

library(ggrepel)

# ── Fig 4: VAF of predicted functional mutations ──
f4 <- ann %>% filter(likely_functional) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol,
                              levels = rev(sort(unique(Hugo_Symbol)))))

p4 <- ggplot(f4, aes(VAF, Hugo_Symbol, color = IMPACT)) +
  geom_vline(xintercept = c(0.1, 0.5), linetype = "dashed", color = "grey60") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = HGVSp_Short), size = 2.6,
                  max.overlaps = Inf, box.padding = 0.4,
                  min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual(values = c(HIGH = "#D7301F", MODERATE = "#2C7BB6"),
                     name = "Impact") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Predicted functional mutations: variant allele fraction",
       subtitle = "Dashed lines at VAF 0.1 and 0.5 (rough guides; VAF depends on tumor purity and copy number)",
       x = "VAF (t_alt_count / t_depth)", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "fig4_functional_VAF.png"),
       p4, width = 9, height = 6, dpi = 300, bg = "white")

print(p3); print(p4)