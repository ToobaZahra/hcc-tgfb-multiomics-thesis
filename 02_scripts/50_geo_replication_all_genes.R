## 50: replication of the 26 nexus genes in GSE14520 and GSE76427
## tests per gene: tumor vs non-tumor, survival (Cox HR per SD), stage (Spearman with TNM; Late vs Early)
## final table puts the TCGA result next to both GEO cohorts. Grade is not available in these cohorts.
library(GEOquery)
library(survival)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
g <- "data/raw/geo/"; a <- "03_results/association"; d_dir <- "03_results/dge"
sig <- function(p) !is.na(p) & p < 0.05

genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
           "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
           "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
           "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# ── helpers ──────────────────────────────────────────────
prep_annot <- function(f) { tab <- Table(getGEO(filename = f))
sc <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
list(ids = as.character(tab$ID), sy = lapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), trimws)) }
gene_vec <- function(ann, mat, gene) {
  idx <- which(vapply(ann$sy, function(x) gene %in% x, logical(1)))
  hit <- intersect(ann$ids[idx], rownames(mat)); if (length(hit) == 0) return(NULL)
  best <- hit[which.max(rowMeans(mat[hit, , drop = FALSE]))]
  list(x = mat[best, ], probe = best, nprobes = length(hit)) }
tnm_num <- function(x) { x <- toupper(trimws(as.character(x)))
ifelse(x == "I", 1, ifelse(x == "II", 2, ifelse(grepl("^III", x), 3, ifelse(grepl("^IV", x), 4, NA)))) }

test_gene <- function(d, mat, ann, gene, ds, linear) {
  gv <- gene_vec(ann, mat, gene)
  if (is.null(gv)) return(data.frame(dataset = ds, gene = gene, note = "no probe"))
  d$x <- unname(gv$x[d$id]); if (linear) d$x <- log2(pmax(d$x, 1))
  r <- data.frame(dataset = ds, gene = gene, probe = gv$probe, n_probes = gv$nprobes, note = "")
  tt <- d[!is.na(d$x) & !is.na(d$tumor), ]
  r$tn_diff <- mean(tt$x[tt$tumor]) - mean(tt$x[!tt$tumor])
  r$tn_p <- wilcox.test(tt$x[tt$tumor], tt$x[!tt$tumor])$p.value
  t <- tt[tt$tumor, ]
  t$stage_grp <- factor(ifelse(t$stage_num %in% 1:2, "Early", ifelse(t$stage_num %in% 3:4, "Late", NA)), c("Early", "Late"))
  s <- t[!is.na(t$time) & t$time > 0 & !is.na(t$event), ]; s$z <- as.numeric(scale(s$x))
  f1 <- tryCatch(summary(coxph(Surv(time, event) ~ z, data = s)), error = function(e) NULL)
  r$HR <- if (is.null(f1)) NA else f1$conf.int[1, 1]
  r$HR_low <- if (is.null(f1)) NA else f1$conf.int[1, 3]; r$HR_high <- if (is.null(f1)) NA else f1$conf.int[1, 4]
  r$cox_p <- if (is.null(f1)) NA else f1$coefficients[1, 5]
  f2 <- tryCatch(summary(coxph(Surv(time, event) ~ z + age + stage_grp, data = s)), error = function(e) NULL)
  r$HR_adj <- if (is.null(f2)) NA else f2$conf.int["z", 1]; r$cox_p_adj <- if (is.null(f2)) NA else f2$coefficients["z", 5]
  k <- t[!is.na(t$stage_num), ]
  ct <- suppressWarnings(cor.test(k$x, k$stage_num, method = "spearman", exact = FALSE))
  r$stage_rho <- unname(ct$estimate); r$stage_p <- ct$p.value
  e <- t[!is.na(t$stage_grp), ]
  r$late_early_diff <- median(e$x[e$stage_grp == "Late"]) - median(e$x[e$stage_grp == "Early"])
  r$late_early_p <- wilcox.test(e$x[e$stage_grp == "Late"], e$x[e$stage_grp == "Early"])$p.value
  r$n_surv <- nrow(s); r$events <- sum(s$event); r
}

# ── cohort 1: GSE14520 ───────────────────────────────────
m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
tis <- tolower(trimws(cl$`Tissue Type`))
d1 <- data.frame(id = as.character(cl$Affy_GSM), tumor = ifelse(grepl("non", tis), FALSE, grepl("tumor", tis)),
                 time = as.numeric(cl$`Survival months`), event = as.numeric(cl$`Survival status`),
                 age = as.numeric(cl$Age), stage_num = tnm_num(cl$`TNM staging`))
d1 <- d1[d1$id %in% colnames(m1), ]
ann1 <- prep_annot(paste0(g, "GPL3921.annot.gz")); lin1 <- max(m1, na.rm = TRUE) > 100

# ── cohort 2: GSE76427 ───────────────────────────────────
m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]]); p2 <- p2[match(colnames(m2), rownames(p2)), ]
tis2 <- tolower(as.character(p2$`tissue:ch1`))
d2 <- data.frame(id = colnames(m2), tumor = ifelse(grepl("non", tis2), FALSE, grepl("tumor", tis2)),
                 time = as.numeric(as.character(p2$`duryears_os:ch1`)) * 12, event = as.numeric(as.character(p2$`event_os:ch1`)),
                 age = as.numeric(as.character(p2$`age (years):ch1`)), stage_num = tnm_num(p2$`tnm_staging_clinical:ch1`))
ann2 <- prep_annot(paste0(g, "GPL10558.annot.gz")); lin2 <- max(m2, na.rm = TRUE) > 100
cat("Linear scale detected: GSE14520 =", lin1, "| GSE76427 =", lin2, "\n")

# ── run all genes ────────────────────────────────────────
res <- bind_rows(lapply(genes, function(gn) bind_rows(
  tryCatch(test_gene(d1, m1, ann1, gn, "GSE14520", lin1), error = function(e) data.frame(dataset = "GSE14520", gene = gn, note = conditionMessage(e))),
  tryCatch(test_gene(d2, m2, ann2, gn, "GSE76427", lin2), error = function(e) data.frame(dataset = "GSE76427", gene = gn, note = conditionMessage(e))))))
for (v in c("tn_p", "cox_p", "stage_p")) res[[paste0("padj_", v)]] <- ave(res[[v]], res$dataset, FUN = function(p) p.adjust(p, "BH"))
write.csv(res, file.path(a, "15_GEO_test_results_long.csv"), row.names = FALSE)
cat("Genes without a probe or with errors:\n"); print(res[res$note != "", c("dataset", "gene", "note")])

# ── TCGA reference and side-by-side table ────────────────
tc_s <- read.csv(file.path(a, "8_expr_vs_survival_gene.csv")) %>% transmute(gene, tcga_HR = HR, tcga_surv_padj = padj_cox)
tc_d <- read.csv(file.path(d_dir, "nexus_stage_grade_dge.csv")) %>%
  transmute(gene, tcga_tn_log2FC = log2FC_tumor_vs_normal, tcga_stage_log2FC = log2FC_stage, tcga_stage_padj = padj_stage, tcga_tn_padj = padj_tumor_vs_normal)
w <- function(ds, suf) res %>% filter(dataset == ds) %>%
  dplyr::select(gene, tn_diff, padj_tn_p, HR, cox_p, padj_cox_p, stage_rho, stage_p, late_early_diff) %>%
  dplyr::rename_with(~ paste0(.x, "_", suf), -gene)
F <- tc_s %>% left_join(tc_d, by = "gene") %>% left_join(w("GSE14520", "A"), by = "gene") %>% left_join(w("GSE76427", "B"), by = "gene")
dirs <- function(x) sign(x)
F$surv_status <- case_when(
  !sig(F$tcga_surv_padj) ~ "TCGA effect not significant (cannot be replicated)",
  dirs(log(F$tcga_HR)) == dirs(log(F$HR_A)) & dirs(log(F$tcga_HR)) == dirs(log(F$HR_B)) & sig(F$cox_p_A) ~ "Replicated (both cohorts same direction, significant in GSE14520)",
  dirs(log(F$tcga_HR)) == dirs(log(F$HR_A)) & dirs(log(F$tcga_HR)) == dirs(log(F$HR_B)) ~ "Same direction in both, not significant",
  TRUE ~ "Not replicated")
F$stage_status <- case_when(
  !sig(F$tcga_stage_padj) ~ "TCGA effect not significant (cannot be replicated)",
  dirs(F$tcga_stage_log2FC) == dirs(F$stage_rho_A) & dirs(F$tcga_stage_log2FC) == dirs(F$stage_rho_B) & sig(F$stage_p_A) ~ "Replicated (both cohorts same direction, significant in GSE14520)",
  dirs(F$tcga_stage_log2FC) == dirs(F$stage_rho_A) & dirs(F$tcga_stage_log2FC) == dirs(F$stage_rho_B) ~ "Same direction in both, not significant",
  TRUE ~ "Not replicated")
F$tn_status <- case_when(
  !sig(F$tcga_tn_padj) ~ "TCGA effect not significant (cannot be replicated)",
  dirs(F$tcga_tn_log2FC) == dirs(F$tn_diff_A) & dirs(F$tcga_tn_log2FC) == dirs(F$tn_diff_B) & sig(F$padj_tn_p_A) & sig(F$padj_tn_p_B) ~ "Replicated in both cohorts",
  dirs(F$tcga_tn_log2FC) == dirs(F$tn_diff_A) & sig(F$padj_tn_p_A) ~ "Replicated in GSE14520 only",
  TRUE ~ "Not replicated")
F <- F %>% arrange(cox_p_A)
write.csv(F, file.path(a, "15_GEO_replication_all_genes.csv"), row.names = FALSE)

cat("\n===== SURVIVAL: TCGA vs GSE14520 (A) vs GSE76427 (B) =====\n")
print(as.data.frame(F %>% dplyr::select(gene, tcga_HR, tcga_surv_padj, HR_A, cox_p_A, HR_B, cox_p_B, surv_status)), digits = 3)
cat("\n===== STAGE =====\n")
print(as.data.frame(F %>% dplyr::select(gene, tcga_stage_log2FC, tcga_stage_padj, stage_rho_A, stage_p_A, stage_rho_B, stage_p_B, stage_status)), digits = 3)
cat("\n===== TUMOR vs NORMAL =====\n")
print(as.data.frame(F %>% dplyr::select(gene, tcga_tn_log2FC, tn_diff_A, padj_tn_p_A, tn_diff_B, padj_tn_p_B, tn_status)), digits = 3)
cat("\nSummary counts:\n"); print(table(F$surv_status)); print(table(F$stage_status)); print(table(F$tn_status))
