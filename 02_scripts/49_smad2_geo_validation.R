## 49: SMAD2 validation in independent expression cohorts (GSE14520, GSE76427)
## Tests: tumor vs non-tumor, survival (Cox + Kaplan-Meier), stage (TNM). Grade is not available in these cohorts.
## TCGA reference: higher expression in tumor (DESeq2 log2FC +0.36); survival HR 1.30 (1.10-1.54) per SD; higher in late stage (DESeq2).
library(GEOquery)
library(survival)
library(dplyr)

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
g <- "data/raw/geo/"; a <- "03_results/association"; fig_dir <- "04_figures"
SYM <- "SMAD2"

# ── helpers ──────────────────────────────────────────────
sym_probe <- function(annot_file, mat, symbol = SYM) {
  tab <- Table(getGEO(filename = annot_file))
  sc  <- grep("symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
  hit <- tab$ID[sapply(strsplit(as.character(tab[[sc]]), " ?/// ?"), function(x) symbol %in% trimws(x))]
  hit <- intersect(as.character(hit), rownames(mat))
  if (length(hit) == 0) stop("No probe for ", symbol, " (symbol column used: ", sc, ")")
  best <- hit[which.max(rowMeans(mat[hit, , drop = FALSE]))]
  cat(annot_file, ": ", length(hit), "probe(s) for", symbol, "-> using", best, "\n"); best
}
tnm_num <- function(x) { x <- toupper(trimws(as.character(x)))
ifelse(x == "I", 1, ifelse(x == "II", 2, ifelse(grepl("^III", x), 3, ifelse(grepl("^IV", x), 4, NA)))) }
out <- data.frame()
add <- function(ds, an, n, eff, p, note = "") out <<- rbind(out, data.frame(dataset = ds, analysis = an, n = n, effect = eff, p = p, note = note))

analyse <- function(ds, d) {          # d: x, tumor (logical), time (months), event, age, stage_num
  cat("\n=====", ds, "=====\n")
  if (any(!is.na(d$tumor))) {
    tt <- d[!is.na(d$x) & !is.na(d$tumor), ]
    add(ds, "Tumor vs non-tumor (mean log2 difference)", nrow(tt),
        mean(tt$x[tt$tumor]) - mean(tt$x[!tt$tumor]), wilcox.test(tt$x[tt$tumor], tt$x[!tt$tumor])$p.value,
        sprintf("%d tumor, %d non-tumor", sum(tt$tumor), sum(!tt$tumor)))
  }
  t <- d[which(d$tumor), ]; t <- t[!is.na(t$x), ]
  t$stage_grp <- factor(ifelse(t$stage_num %in% 1:2, "Early", ifelse(t$stage_num %in% 3:4, "Late", NA)), c("Early", "Late"))
  cat("Tumors:", nrow(t), "| stage groups:\n"); print(table(t$stage_grp, useNA = "ifany"))
  # survival
  s <- t[!is.na(t$time) & t$time > 0 & !is.na(t$event), ]; s$z <- as.numeric(scale(s$x))
  f1 <- coxph(Surv(time, event) ~ z, data = s); u <- summary(f1)
  add(ds, "Survival: Cox HR per 1 SD", f1$n, u$conf.int[1, 1], u$coefficients[1, 5],
      sprintf("95%% CI %.2f-%.2f; %d events", u$conf.int[1, 3], u$conf.int[1, 4], f1$nevent))
  f2 <- tryCatch(coxph(Surv(time, event) ~ z + age + stage_grp, data = s), error = function(e) NULL)
  if (!is.null(f2)) { v <- summary(f2)
  add(ds, "Survival: Cox HR adjusted for age + stage", f2$n, v$conf.int["z", 1], v$coefficients["z", 5],
      sprintf("95%% CI %.2f-%.2f; %d events", v$conf.int["z", 3], v$conf.int["z", 4], f2$nevent)) }
  s$grp <- factor(ifelse(s$x > median(s$x), "High", "Low"), c("Low", "High"))
  lr <- survdiff(Surv(time, event) ~ grp, data = s)
  add(ds, "Survival: median split log-rank", nrow(s), NA, 1 - pchisq(lr$chisq, 1), "")
  plot(survfit(Surv(time, event) ~ grp, data = s), col = c("#2C7BB6", "#D7301F"), lwd = 2, xlab = "Months", ylab = "Overall survival",
       main = sprintf("%s: %s\nlog-rank p = %.3g", ds, SYM, 1 - pchisq(lr$chisq, 1)), cex.main = 0.9)
  legend("bottomleft", c("Low", "High"), col = c("#2C7BB6", "#D7301F"), lwd = 2, bty = "n", cex = 0.8)
  # stage
  k <- t[!is.na(t$stage_num), ]
  if (nrow(k) > 20) { ct <- suppressWarnings(cor.test(k$x, k$stage_num, method = "spearman", exact = FALSE))
  add(ds, "Stage (TNM I-IV): Spearman rho", nrow(k), unname(ct$estimate), ct$p.value, paste(names(table(k$stage_num)), table(k$stage_num), sep = "=", collapse = " ")) }
  e <- t[!is.na(t$stage_grp), ]
  if (nrow(e) > 20) add(ds, "Late vs Early stage (median log2 difference)", nrow(e),
                        median(e$x[e$stage_grp == "Late"]) - median(e$x[e$stage_grp == "Early"]),
                        wilcox.test(e$x[e$stage_grp == "Late"], e$x[e$stage_grp == "Early"])$p.value,
                        sprintf("%d early, %d late", sum(e$stage_grp == "Early"), sum(e$stage_grp == "Late")))
}

png(file.path(fig_dir, "fig_SMAD2_GEO_km.png"), width = 2000, height = 1000, res = 250)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# ── GSE14520 (Affymetrix GPL3921; clinical in the supplementary file) ──
m1 <- readRDS(paste0(g, "GSE14520_norm.rds"))
cl <- read.delim(gzfile(paste0(g, "GSE14520_Extra_Supplement.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
pr1 <- sym_probe(paste0(g, "GPL3921.annot.gz"), m1)
x1 <- m1[pr1, ]
tis <- tolower(trimws(cl$`Tissue Type`)); cat("Tissue types:\n"); print(table(cl$`Tissue Type`))
cat("TNM staging values:\n"); print(table(cl$`TNM staging`, useNA = "ifany"))
d1 <- data.frame(x = unname(x1[as.character(cl$Affy_GSM)]), tumor = ifelse(grepl("non", tis), FALSE, grepl("tumor", tis)),
                 time = as.numeric(cl$`Survival months`), event = as.numeric(cl$`Survival status`),
                 age = as.numeric(cl$Age), stage_num = tnm_num(cl$`TNM staging`))
if (max(m1, na.rm = TRUE) > 100) { cat("GSE14520 matrix is on a linear scale: applying log2 to the SMAD2 probe\n"); d1$x <- log2(pmax(d1$x, 1)) }
analyse("GSE14520", d1)

# ── GSE76427 (Illumina GPL10558; clinical in the sample table) ──
m2 <- readRDS(paste0(g, "GSE76427_norm.rds"))
p2 <- Biobase::pData(readRDS(paste0(g, "GSE76427.rds"))[[1]])
p2 <- p2[match(colnames(m2), rownames(p2)), ]
pr2 <- sym_probe(paste0(g, "GPL10558.annot.gz"), m2)
tis2 <- tolower(as.character(p2$`tissue:ch1`)); cat("Tissue values:\n"); print(table(p2$`tissue:ch1`))
d2 <- data.frame(x = unname(m2[pr2, ]), tumor = ifelse(grepl("non", tis2), FALSE, grepl("tumor", tis2)),
                 time = as.numeric(as.character(p2$`duryears_os:ch1`)) * 12, event = as.numeric(as.character(p2$`event_os:ch1`)),
                 age = as.numeric(as.character(p2$`age (years):ch1`)), stage_num = tnm_num(p2$`tnm_staging_clinical:ch1`))
if (max(m2, na.rm = TRUE) > 100) { cat("GSE76427 matrix is on a linear scale: applying log2 to the SMAD2 probe\n"); d2$x <- log2(pmax(d2$x, 1)) }
analyse("GSE76427", d2)
dev.off()

out$p <- signif(out$p, 3); out$effect <- signif(out$effect, 3)
write.csv(out, file.path(a, "14_SMAD2_GEO_validation.csv"), row.names = FALSE)
cat("\n===== SMAD2 validation summary =====\n"); print(out, row.names = FALSE, right = FALSE)
cat("\nTCGA-LIHC reference: Cox HR 1.30 (1.10-1.54), p = 0.0026; DESeq2 log2FC tumor vs normal +0.36; Late vs Early log2FC +0.16 (padj 0.033)\n")