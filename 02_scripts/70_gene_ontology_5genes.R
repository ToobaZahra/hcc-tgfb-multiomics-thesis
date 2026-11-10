## 70: Gene Ontology (GO) annotations of DUSP10, SMAD2, MYC, MAPK8, HDAC11
## Source: org.Hs.eg.db (GO annotations of each gene) and GO.db (term names). No internet needed.
## Outputs (03_results/pathways):
##   40_GO_annotations_5genes.csv       every GO term of each gene (BP, MF, CC) with the evidence code
##   40b_GO_themes_by_gene.csv          counts of GO terms per theme and gene (readable summary)
##   40c_GO_shared_terms.csv            GO terms shared by 2 or more of the genes
##   fig_5genes_GO_themes.png           tile plot: theme x gene
## Note: these are the terms DIRECTLY annotated to each gene (parent terms are not added).
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(GO.db); library(dplyr); library(tidyr); library(ggplot2) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
out_dir <- "03_results/pathways"; fig_dir <- "04_figures"
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")          # add "NPC1" here if you want it too

# ---- 1. GO annotations of the genes ----
ann <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
                             columns = c("ENTREZID", "GENENAME", "GO", "EVIDENCE", "ONTOLOGY"))
ann <- ann[!is.na(ann$GO), ]
cat("Gene names:\n"); print(unique(ann[, c("SYMBOL", "ENTREZID", "GENENAME")]), row.names = FALSE)

terms <- AnnotationDbi::select(GO.db, keys = unique(ann$GO), keytype = "GOID", columns = c("TERM"))
ann <- left_join(ann, terms, by = c("GO" = "GOID"))
exp_codes <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")
ann$evidence_class <- ifelse(ann$EVIDENCE %in% exp_codes, "experimental",
                             ifelse(ann$EVIDENCE == "IEA", "electronic (IEA)", "curated / inferred"))
ann$ONTOLOGY <- factor(ann$ONTOLOGY, c("BP", "MF", "CC"), c("Biological process", "Molecular function", "Cellular component"))
ann <- ann %>% dplyr::select(gene = SYMBOL, ontology = ONTOLOGY, go_id = GO, term = TERM, evidence = EVIDENCE, evidence_class) %>%
  distinct() %>% arrange(gene, ontology, term)
write.csv(ann, file.path(out_dir, "40_GO_annotations_5genes.csv"), row.names = FALSE)

cat("\nNumber of GO terms per gene and ontology (all evidence, and experimental / curated only):\n")
print(ann %>% group_by(gene, ontology) %>% summarise(all = n_distinct(go_id),
                                                     non_IEA = n_distinct(go_id[evidence_class != "electronic (IEA)"]), .groups = "drop") %>%
        pivot_wider(names_from = ontology, values_from = c(all, non_IEA)), width = Inf)

# ---- 2. Terms shared by two or more genes (non-IEA, so the list is more reliable) ----
core <- ann %>% filter(evidence_class != "electronic (IEA)")
shared <- core %>% group_by(ontology, go_id, term) %>% summarise(n_genes = n_distinct(gene),
                                                                 genes = paste(sort(unique(gene)), collapse = ", "), .groups = "drop") %>%
  filter(n_genes >= 2) %>% arrange(desc(n_genes), ontology, term)
write.csv(shared, file.path(out_dir, "40c_GO_shared_terms.csv"), row.names = FALSE)
cat("\n== GO terms shared by 2 or more genes (non-IEA) ==\n"); print(as.data.frame(shared), right = FALSE, row.names = FALSE)

# ---- 3. Readable themes (keyword groups over the term names) ----
themes <- c(
  "MAPK / JNK / p38 cascade"       = "map kinase|mapk|jnk|p38|stress-activated",
  "TGF-beta / SMAD / BMP signalling" = "transforming growth factor|tgf|smad|activin|nodal|bmp",
  "Cell cycle and proliferation"   = "cell cycle|proliferation|mitotic|cell division|cell growth",
  "Apoptosis / cell death"         = "apopto|cell death|anoikis",
  "Regulation of transcription"    = "transcription|rna polymerase ii",
  "Chromatin / histone"            = "histone|chromatin|deacetylat|acetylat",
  "Phosphatase / kinase / phosphorylation" = "phosphatase|kinase|dephosphorylation|phosphorylation",
  "Stress, inflammation and immunity" = "stress|inflammat|immune|cytokine|interleukin|lipopolysaccharide|toll-like",
  "Development and differentiation" = "development|differentiation|morphogenesis")
tl <- bind_rows(lapply(names(themes), function(th)
  core %>% filter(grepl(themes[[th]], term, ignore.case = TRUE)) %>% mutate(theme = th)))
th_count <- tl %>% group_by(theme, gene) %>% summarise(n_terms = n_distinct(go_id), .groups = "drop") %>%
  complete(theme = names(themes), gene = genes, fill = list(n_terms = 0))
write.csv(th_count %>% pivot_wider(names_from = gene, values_from = n_terms), file.path(out_dir, "40b_GO_themes_by_gene.csv"), row.names = FALSE)
cat("\n== GO terms per theme and gene (non-IEA) ==\n")
print(as.data.frame(th_count %>% pivot_wider(names_from = gene, values_from = n_terms)), right = FALSE, row.names = FALSE)

th_count$theme <- factor(th_count$theme, rev(names(themes))); th_count$gene <- factor(th_count$gene, genes)
p <- ggplot(th_count, aes(gene, theme, fill = n_terms)) + geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n_terms > 0, n_terms, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "#2C7BB6", name = "GO terms") +
  labs(title = "Gene Ontology themes of the five genes",
       subtitle = "Number of directly annotated GO terms (experimental / curated evidence) matching each theme", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "fig_5genes_GO_themes.png"), p, width = 8, height = 5, dpi = 300, bg = "white")
cat("\nSaved 40_, 40b_, 40c_ tables to", out_dir, "and fig_5genes_GO_themes.png to", fig_dir, "\n")
cat("Themes are keyword groups over GO term names, so they are a reading aid, not an enrichment test.\n")