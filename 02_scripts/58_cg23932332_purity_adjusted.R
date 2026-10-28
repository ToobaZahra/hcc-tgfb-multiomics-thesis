## 58: Is the cg23932332 vs DUSP10 expression correlation (rho = -0.40) explained by tumor purity?
## Normal liver cells have lower methylation at this CpG and higher DUSP10 expression, so tumors with more
## normal-cell contamination could show a negative correlation without any silencing.
## Tests: raw Spearman, purity-adjusted partial Spearman (several purity estimates), and by grade group.
suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(org.Hs.eg.db) })

setwd("E:/MS-Project/hcc-tgfb-multiomics-thesis")
a <- "03_results/association"
G <- "DUSP10"; key <- "cg23932332"

# ---- 1. TCGA-LIHC tumors: probe beta, DUSP10 expression, grade ----
se <- readRDS("data/raw/tcga_lihc_meth450_se.rds")
bc <- substr(colnames(se), 14, 15)
tum <- which(bc == "01"); tum <- tum[!duplicated(substr(colnames(se)[tum], 1, 12))]
cdm <- as.data.frame(colData(se))[tum, ]
d <- data.frame(patient = substr(colnames(se)[tum], 1, 12),
                sample15 = substr(colnames(se)[tum], 1, 15),
                beta = as.numeric(assay(se[key, tum])),
                grade = factor(cdm$tumor_grade, c("G1", "G2", "G3", "G4")))

vsd <- readRDS("data/raw/tcga_lihc_vsd.rds"); expr <- assay(vsd)
ens <- suppressMessages(mapIds(org.Hs.eg.db, G, "ENSEMBL", "SYMBOL", multiVals = "first"))
row <- match(ens, sub("\\..*$", "", rownames(expr))); stopifnot(!is.na(row))
et <- which(substr(colnames(expr), 14, 15) == "01"); et <- et[!duplicated(substr(colnames(expr)[et], 1, 12))]
ex <- expr[row, et]; names(ex) <- substr(colnames(expr)[et], 1, 12)
d$expr <- unname(ex[d$patient])

# ---- 2. Tumor purity estimates (TCGAbiolinks ships the Aran et al. 2015 consensus estimates) ----
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) stop("Install TCGAbiolinks (BiocManager::install('TCGAbiolinks')) to get purity estimates.")
pur <- TCGAbiolinks::Tumor.purity
cat("Tumor.purity columns:", paste(names(pur), collapse = ", "), "\n"); print(head(pur, 3))
idcol <- grep("sample|barcode|^id", names(pur), ignore.case = TRUE, value = TRUE)[1]
cat("Using ID column:", idcol, "| example ID:", as.character(pur[[idcol]][1]), "| example from our data:", d$sample15[1], "\n")
ids <- as.character(pur[[idcol]])
pur$patient <- substr(ids, 1, 12)                           # match on patient, not on the full sample barcode
keep <- substr(ids, 14, 15) %in% "01"; if (sum(keep) > 0) pur <- pur[keep, ]   # primary tumor rows only
num <- function(x) as.numeric(sub(",", ".", as.character(x)))   # some versions store decimals with commas
pcols <- intersect(c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"), names(pur))
for (cc in pcols) pur[[cc]] <- num(pur[[cc]])
pur <- pur[!duplicated(pur$patient), c("patient", pcols)]
d <- left_join(d, pur, by = "patient")
metrics <- intersect(c("CPE", "ABSOLUTE", "ESTIMATE", "LUMP", "IHC"), names(d))
cat("Tumors:", nrow(d), "| with purity per metric:\n"); print(sapply(d[metrics], function(x) sum(is.finite(x))))
if (all(sapply(d[metrics], function(x) sum(is.finite(x))) == 0))
  stop("No purity values matched. Paste the 'Tumor.purity columns' line and the example IDs printed above.")
metrics <- metrics[sapply(d[metrics], function(x) sum(is.finite(x))) > 0]

# ---- 3. Tests ----
sp <- function(x, y) { k <- is.finite(x) & is.finite(y)
ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE))
c(rho = unname(ct$estimate), p = ct$p.value, n = sum(k)) }
# partial Spearman: rank-transform, remove the linear effect of rank(purity), correlate the residuals
psp <- function(x, y, z) { k <- is.finite(x) & is.finite(y) & is.finite(z)
rx <- resid(lm(rank(x[k]) ~ rank(z[k]))); ry <- resid(lm(rank(y[k]) ~ rank(z[k])))
r <- cor(rx, ry); df <- sum(k) - 3; t <- r * sqrt(df / (1 - r^2))
c(rho = r, p = 2 * pt(-abs(t), df), n = sum(k)) }
row1 <- function(label, v) data.frame(analysis = label, rho = v["rho"], p = v["p"], n = v["n"], row.names = NULL)

res <- bind_rows(
  row1("Raw: beta vs expression (all tumors)", sp(d$beta, d$expr)),
  bind_rows(lapply(metrics, function(m) row1(paste0("Adjusted for purity (", m, ")"), psp(d$beta, d$expr, d[[m]])))),
  bind_rows(lapply(metrics, function(m) row1(paste0("Purity (", m, ") vs beta"), sp(d[[m]], d$beta)))),
  bind_rows(lapply(metrics, function(m) row1(paste0("Purity (", m, ") vs expression"), sp(d[[m]], d$expr)))))

# by grade group (raw)
d$grade_grp <- factor(ifelse(d$grade %in% c("G1", "G2"), "G1-G2", ifelse(d$grade %in% c("G3", "G4"), "G3-G4", NA)))
strat <- bind_rows(lapply(levels(d$grade_grp), function(g) {
  s <- d[which(d$grade_grp == g), ]; row1(paste0("Raw within ", g), sp(s$beta, s$expr)) }))
res <- bind_rows(res, strat)

write.csv(res, file.path(a, "28_cg23932332_purity_adjusted.csv"), row.names = FALSE)
cat("\n== cg23932332 vs DUSP10 expression, purity check ==\n"); print(res, digits = 3)
cat("\nHow to read it: if the adjusted rho stays near -0.40 (and p stays small), purity does not explain the link.\n")
cat("If it falls toward 0, the correlation was largely tumor/normal mixing. Purity vs beta and vs expression show the two legs of that confounding.\n")
cat("The grade strata check whether the link holds within comparable tumors.\n")