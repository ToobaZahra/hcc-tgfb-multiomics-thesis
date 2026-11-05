## 68: Stage, grade, age and sex versus expression and methylation for DUSP10, SMAD2, MYC, MAPK8, HDAC11
##   Expression: TCGA-LIHC (stage, grade, age, sex), GSE14520 and GSE76427 (stage, age; no grade in these files)
##   Methylation (gene-level mean beta and strongest probe): TCGA-LIHC only (the other methylation cohorts have no clinical data)
## Tests: Spearman rho for the ordered features (stage I-IV, grade G1-G4, age); Wilcoxon for sex. Raw p-values.
## Figures go to 04_figures; every number goes to 03_results/association/38_clinical_associations_5genes.csv
suppressPackageStartupMessages({ library(SummarizedExperiment); library(DESeq2); library(org.Hs.eg.db)
  library(GEOquery); library(ggplot2); library(dplyr) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
fig_dir <- "04_figures"; a <- "03_results/association"; g <- "data/raw/geo/"
dir.create(fig_dir, showWarnings = FALSE)
genes <- c("DUSP10", "SMAD2", "MYC", "MAPK8", "HDAC11")

# ================= 1. TCGA clinical features (same fields as scripts 51 and 62) =================
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cd <- as.data.frame(colData(se))[tum, ]
patient <- substr(colnames(se)[tum], 1, 12)
stx <- as.character(cd$ajcc_pathologic_stage)
stage_num <- ifelse(grepl("^Stage IV", stx), 4, ifelse(grepl("^Stage III", stx), 3,
                                                       ifelse(grepl("^Stage II", stx), 2, ifelse(grepl("^Stage I", stx), 1, NA))))
grade_num <- as.integer(factor(cd$tumor_grade, c("G1", "G2", "G3", "G4")))
agecol <- intersect(c("age_at_index", "age_at_diagnosis"), names(cd))[1]
age <- as.numeric(cd[[agecol]]); if (median(age, na.rm = TRUE) > 200) age <- age / 365.25
sexcol <- intersect(c("gender", "sex"), names(cd))[1]
sex <- if (is.na(sexcol)) rep(NA_character_, nrow(cd)) else tolower(as.character(cd[[sexcol]]))
clin_t <- data.frame(patient = patient, stage_num = stage_num, grade_num = grade_num, age = age, sex = sex)
cat("TCGA tumors:", nrow(clin_t), "| stage known:", sum(!is.na(stage_num)), "| grade known:", sum(!is.na(grade_num)),
    "| age known:", sum(is.finite(age)), "| sex values:", paste(names(table(sex)), collapse = "/"), "\n")

# ================= 2. TCGA expression and methylation per gene =================
vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
pat_e <- substr(colnames(expr)[et], 1, 12)
ens <- suppressMessages(mapIds(org.Hs.eg.db, genes, "ENSEMBL", "SYMBOL", multiVals = "first"))
rw <- match(ens, sub("\\..*$", "", rownames(expr)))
if (any(is.na(rw))) stop("Not found in the TCGA expression matrix: ", paste(genes[is.na(rw)], collapse = ", "))

pr_all <- read.csv("03_results/methylation/meth_expr_corr_probe_level.csv")
stopifnot(all(c("gene", "probe", "rho") %in% names(pr_all)))
keyp <- character(0); gprobes <- list()
for (G in genes) {
  pg <- pr_all[pr_all$gene == G & pr_all$probe %in% rownames(se), ]
  if (nrow(pg) == 0) { cat("No methylation probes for", G, ": skipped for methylation.\n"); next }
  gprobes[[G]] <- pg$probe
  keyp[G] <- if (G == "DUSP10" && "cg23932332" %in% pg$probe) "cg23932332" else pg$probe[which.max(abs(pg$rho))]
}
mg <- names(gprobes)

mk_long <- function(cohort, gene, measure, ids, vals, clin) {
  d <- data.frame(patient = ids, value = as.numeric(vals))
  d <- inner_join(d, clin, by = "patient"); d$cohort <- cohort; d$gene <- gene; d$measure <- measure; d }
tc_expr <- bind_rows(lapply(seq_along(genes), function(i) mk_long("TCGA-LIHC", genes[i], "Expression", pat_e, expr[rw[i], et], clin_t)))
tc_mg <- bind_rows(lapply(mg, function(G) mk_long("TCGA-LIHC", G, "Gene-level methylation", patient,
                                                  colMeans(assay(se[gprobes[[G]], tum]), na.rm = TRUE), clin_t)))
tc_mk <- bind_rows(lapply(mg, function(G) mk_long("TCGA-LIHC", paste0(G, "\n", keyp[[G]]), "Strongest-probe methylation", patient,
                                                  assay(se[keyp[[G]], tum]), clin_t)))

# ================= 3. GEO expression cohorts: stage and age (tumors only) =================
prep_annot <- function(f) { tab <- Table(getGEO(filename = f))
sc <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
list(ids = as.character(tab$ID), sy = lapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), trimws)) }
gene_vec <- function(ann, mat, gene) {
  idx <- which(vapply(ann$sy, function(x) gene %in% x, logical(1)))
  hit <- intersect(ann$ids[idx], rownames(mat)); if (length(hit) == 0) return(NULL)
  best <- hit[which.max(rowMeans(mat[hit, , drop = FALSE]))]; mat[best, ] }
tnm_num <- function(x) { x <- toupper(trimws(as.character(x)))
ifelse(x == "I", 1, ifelse(x == "II", 2, ifelse(grepl("^III", x), 3, ifelse(grepl("^IV", x), 4, NA)))) }

m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
tis <- tolower(trimws(cl$`Tissue Type`))
c1 <- data.frame(id = as.character(cl$Affy_GSM), tumor = !grepl("non", tis) & grepl("tumor", tis),
                 stage_num = tnm_num(cl$`TNM staging`), age = as.numeric(cl$Age))
c1 <- c1[c1$tumor & c1$id %in% colnames(m1) & !duplicated(c1$id), ]
ann1 <- prep_annot(paste0(g, "GPL3921.annot.gz")); lin1 <- max(m1, na.rm = TRUE) > 100
m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]]); p2 <- p2[match(colnames(m2), rownames(p2)), ]
tis2 <- tolower(as.character(p2$`tissue:ch1`))
c2 <- data.frame(id = colnames(m2), tumor = !grepl("non", tis2) & grepl("tumor", tis2),
                 stage_num = tnm_num(p2$`tnm_staging_clinical:ch1`), age = as.numeric(as.character(p2$`age (years):ch1`)))
c2 <- c2[c2$tumor & !duplicated(c2$id), ]
ann2 <- prep_annot(paste0(g, "GPL10558.annot.gz")); lin2 <- max(m2, na.rm = TRUE) > 100
cat("GSE14520 tumors:", nrow(c1), "| GSE76427 tumors:", nrow(c2), "\n")

geo_long <- function(cohort, cdf, m, ann, lin, G) {
  v <- gene_vec(ann, m, G); if (is.null(v)) { cat(G, "not found on", cohort, "\n"); return(NULL) }
  x <- unname(v[cdf$id]); if (lin) x <- log2(pmax(x, 1))
  data.frame(patient = cdf$id, value = x, stage_num = cdf$stage_num, grade_num = NA_integer_, age = cdf$age, sex = NA_character_,
             cohort = cohort, gene = G, measure = "Expression") }
geo_expr <- bind_rows(lapply(genes, function(G) bind_rows(geo_long("GSE14520", c1, m1, ann1, lin1, G),
                                                          geo_long("GSE76427", c2, m2, ann2, lin2, G))))
expr_all <- bind_rows(tc_expr, geo_expr)

# ================= 4. Statistics =================
sp <- function(x, y) { k <- is.finite(x) & is.finite(y)
if (sum(k) < 10 || length(unique(y[k])) < 2) return(c(est = NA, p = NA, n = sum(k)))
ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE))
c(est = unname(ct$estimate), p = ct$p.value, n = sum(k)) }
sex_test <- function(x, s) { k <- is.finite(x) & s %in% c("male", "female")
if (sum(k & s == "male") < 5 || sum(k & s == "female") < 5) return(c(est = NA, p = NA, n = sum(k)))
c(est = median(x[k & s == "male"]) - median(x[k & s == "female"]),
  p = wilcox.test(x[k & s == "male"], x[k & s == "female"])$p.value, n = sum(k)) }
all_rows <- bind_rows(expr_all, tc_mg, tc_mk)
stats <- bind_rows(lapply(split(all_rows, list(all_rows$cohort, all_rows$gene, all_rows$measure), drop = TRUE), function(d) {
  mk <- function(feat, stat, v) data.frame(cohort = d$cohort[1], gene = sub("\n.*", "", d$gene[1]), probe = ifelse(grepl("\n", d$gene[1]), sub(".*\n", "", d$gene[1]), NA),
                                           measure = d$measure[1], feature = feat, statistic = stat, estimate = unname(v["est"]), p = unname(v["p"]), n = unname(v["n"]))
  rbind(mk("Stage I-IV", "Spearman rho", sp(d$value, d$stage_num)), mk("Grade G1-G4", "Spearman rho", sp(d$value, d$grade_num)),
        mk("Age", "Spearman rho", sp(d$value, d$age)), mk("Sex (male minus female)", "median difference", sex_test(d$value, d$sex))) }))
stats <- stats[!is.na(stats$estimate), ]
write.csv(stats, file.path(a, "38_clinical_associations_5genes.csv"), row.names = FALSE)

# ================= 5. Figures =================
lab_of <- function(feature, by) { s <- stats[stats$feature == feature, ]
s$gene_row <- if (by == "probe") paste0(s$gene, "\n", s$probe) else s$gene
s$label <- sprintf("rho = %.2f, p = %.2g\nn = %d", s$estimate, s$p, s$n); s }
pal <- scale_fill_brewer(palette = "Blues", guide = "none")
sub_txt <- "Spearman rho (raw p) for the ordered feature; boxes = median and quartiles"

box_fig <- function(df, xvar, xlevels, xlabels, row_col, col_col, rowlev, collev, feat, by, ylab, title, file, w, h) {
  df$x <- factor(df[[xvar]], xlevels, xlabels); df <- df[!is.na(df$x) & is.finite(df$value), ]
  df$row <- factor(df[[row_col]], rowlev); df$col <- factor(df[[col_col]], collev)
  lb <- lab_of(feat, by); lb$row <- factor(lb[[ifelse(by == "probe", "gene_row", "gene")]], rowlev); lb$col <- factor(lb$cohort, collev)
  lb <- lb[lb$measure %in% unique(df$measure) & !is.na(lb$row) & !is.na(lb$col), ]
  p <- ggplot(df, aes(x, value, fill = x)) + geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey20") +
    geom_text(data = lb, aes(x = 1, y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 2.7) +
    facet_grid(row ~ col, scales = "free_y") + pal + scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(title = title, subtitle = sub_txt, x = NULL, y = ylab) +
    theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, file), p, width = w, height = h, dpi = 300, bg = "white") }

cohE <- c("TCGA-LIHC", "GSE14520", "GSE76427")
ee <- expr_all; ee$stage_f <- ee$stage_num; ee$grade_f <- ee$grade_num
box_fig(ee, "stage_f", 1:4, c("I", "II", "III", "IV"), "gene", "cohort", genes, cohE, "Stage I-IV", "gene",
        "Expression (TCGA: VST; GEO: log2)", "Expression by stage", "fig_5genes_expression_by_stage.png", 10, 12)
box_fig(tc_expr %>% mutate(grade_f = grade_num), "grade_f", 1:4, c("G1", "G2", "G3", "G4"), "gene", "cohort", genes, "TCGA-LIHC", "Grade G1-G4", "gene",
        "Expression (VST)", "Expression by grade (TCGA-LIHC)", "fig_5genes_expression_by_grade.png", 4.2, 12)

lv_stage <- rbind(tc_mg %>% mutate(g2 = gene), tc_mk %>% mutate(g2 = gene)); lv_stage$stage_f <- lv_stage$stage_num; lv_stage$grade_f <- lv_stage$grade_num
box_fig(tc_mg %>% mutate(stage_f = stage_num), "stage_f", 1:4, c("I", "II", "III", "IV"), "gene", "cohort", mg, "TCGA-LIHC", "Stage I-IV", "gene",
        "Mean beta (gene's probes)", "Gene-level methylation by stage (TCGA-LIHC)", "fig_5genes_methylation_by_stage.png", 4.2, 12)
box_fig(tc_mg %>% mutate(grade_f = grade_num), "grade_f", 1:4, c("G1", "G2", "G3", "G4"), "gene", "cohort", mg, "TCGA-LIHC", "Grade G1-G4", "gene",
        "Mean beta (gene's probes)", "Gene-level methylation by grade (TCGA-LIHC)", "fig_5genes_methylation_by_grade.png", 4.2, 12)
pl_rows <- paste0(mg, "\n", keyp[mg])
box_fig(tc_mk %>% mutate(stage_f = stage_num), "stage_f", 1:4, c("I", "II", "III", "IV"), "gene", "cohort", pl_rows, "TCGA-LIHC", "Stage I-IV", "probe",
        "Beta (one probe per gene)", "Strongest-probe methylation by stage (TCGA-LIHC)", "fig_5genes_methylation_key_probe_by_stage.png", 4.2, 12)
box_fig(tc_mk %>% mutate(grade_f = grade_num), "grade_f", 1:4, c("G1", "G2", "G3", "G4"), "gene", "cohort", pl_rows, "TCGA-LIHC", "Grade G1-G4", "probe",
        "Beta (one probe per gene)", "Strongest-probe methylation by grade (TCGA-LIHC)", "fig_5genes_methylation_key_probe_by_grade.png", 4.2, 12)

# age: scatter with a smooth line
age_fig <- function(df, rowlev, collev, by, ylab, title, file, w, h) {
  df <- df[is.finite(df$age) & is.finite(df$value), ]; df$row <- factor(df$gene, rowlev); df$col <- factor(df$cohort, collev)
  lb <- lab_of("Age", by); lb$row <- factor(lb[[ifelse(by == "probe", "gene_row", "gene")]], rowlev); lb$col <- factor(lb$cohort, collev)
  lb <- lb[lb$measure %in% unique(df$measure) & !is.na(lb$row) & !is.na(lb$col), ]
  p <- ggplot(df, aes(age, value)) + geom_point(size = 0.6, alpha = 0.35, color = "grey30") +
    geom_smooth(method = "lm", se = FALSE, color = "#D7301F", linewidth = 0.6) +
    geom_text(data = lb, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.15, size = 2.7) +
    facet_grid(row ~ col, scales = "free_y") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(title = title, subtitle = "Spearman rho (raw p); red line = linear fit", x = "Age (years)", y = ylab) +
    theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, file), p, width = w, height = h, dpi = 300, bg = "white") }
age_fig(expr_all, genes, cohE, "gene", "Expression (TCGA: VST; GEO: log2)", "Expression by age", "fig_5genes_expression_by_age.png", 10, 12)
age_fig(tc_mg, mg, "TCGA-LIHC", "gene", "Mean beta (gene's probes)", "Gene-level methylation by age (TCGA-LIHC)", "fig_5genes_methylation_by_age.png", 4.2, 12)
age_fig(tc_mk, pl_rows, "TCGA-LIHC", "probe", "Beta (one probe per gene)", "Strongest-probe methylation by age (TCGA-LIHC)", "fig_5genes_methylation_key_probe_by_age.png", 4.2, 12)

# sex (TCGA only): expression and methylation, genes across the columns
sx <- bind_rows(tc_expr, tc_mg, tc_mk %>% mutate(gene = sub("\n.*", "", gene)))
sx <- sx[sx$sex %in% c("male", "female") & is.finite(sx$value), ]
sx$measure <- factor(sx$measure, c("Expression", "Gene-level methylation", "Strongest-probe methylation")); sx$gene <- factor(sx$gene, genes)
slab <- stats[stats$feature == "Sex (male minus female)", ]; slab$measure <- factor(slab$measure, levels(sx$measure)); slab$gene <- factor(slab$gene, genes)
slab$label <- sprintf("diff = %.2g, p = %.2g", slab$estimate, slab$p)
ps <- ggplot(sx, aes(sex, value, fill = sex)) + geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, color = "grey20") +
  geom_text(data = slab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.15, size = 2.6) +
  facet_grid(measure ~ gene, scales = "free_y") + scale_fill_manual(values = c(female = "#C51B8A", male = "#2C7BB6"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(title = "Expression and methylation by sex (TCGA-LIHC)", subtitle = "Wilcoxon test (raw p); difference = male minus female median", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_5genes_by_sex_TCGA.png"), ps, width = 14, height = 8, dpi = 300, bg = "white")

cat("\n== Clinical associations (raw p-values) ==\n"); print(as.data.frame(stats), digits = 3)
cat("\nSaved figures to", fig_dir, "and 38_clinical_associations_5genes.csv\n")
cat("Methylation cohorts GSE56588 and GSE54503 have no clinical data, so methylation vs stage/grade/age/sex is TCGA only.\n")
cat("GEO files used here have stage and age but no grade. The strongest probe was chosen on TCGA expression correlation.\n")