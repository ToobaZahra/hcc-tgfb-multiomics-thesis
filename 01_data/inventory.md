# TCGA-LIHC + GEO Data Inventory
**Thesis:** Multi-Omics Analysis of TGF-β Pathway Dysregulation in Hepatocellular Carcinoma
**Author:** Tooba Zahra | MS Bioinformatics, COMSATS University Islamabad
**Last updated:** 24 April 2026

---

## TCGA-LIHC (Primary Cohort)
**Source:** GDC Portal (verified via TCGAbiolinks R package)
**Role:** Primary analysis — DGE, methylation, mutation, staging, survival

| Data Type | n samples | Notes |
|---|---|---|
| Primary Tumor (RNA-seq, STAR-Counts) | 371 | TCGAbiolinks query verified 28 Apr 2026 |
| Solid Tissue Normal (RNA-seq, STAR-Counts) | 50 | TCGAbiolinks query verified 28 Apr 2026 |
| DNA Methylation (450K) | 379 | Illumina 450K arrays |
| Somatic Mutation (MAF) | 373 | Masked somatic mutations |
| Clinical (staging) | 377 | Full AJCC staging available |
| Proteome Profiling (RPPA) | 184 | Bonus — check for phospho-SMAD antibodies |

**Total RNA-seq for DESeq2:** 371 tumor + 50 normal = **421 samples**

---

## GEO Replication Cohorts

### GSE14520
**Role:** DGE replication + stage-stratified replication
**Platform:** Affymetrix HG-U133A 2.0 (GPL3921 primary, GPL571 secondary)
**Source:** Roessler et al., 2010

| Sample Type | n |
|---|---|
| Tumor | 247 |
| Non-tumor | 241 |
| **Total** | **488** |

- Staging: ✅ **TNM staging available** (in clinical supplementary file)
- Note: Dual platform — use GPL3921 (n=445) as primary, GPL571 (n=43) as secondary
- Requires batch correction if both platforms combined

---

### GSE36376
**Role:** DGE replication only (no staging)
**Platform:** Illumina HumanHT-12 V4.0 (GPL10558)
**Source:** Lim et al., 2013

| Sample Type | n |
|---|---|
| Tumor | 240 |
| Non-tumor | 193 |
| **Total** | **433** |

- Staging: ❌ Not available in GEO supplementary (only raw + non-normalized expression files)
- Use only for tumor vs normal DGE replication

---

### GSE76427
**Role:** DGE replication + stage-stratified replication
**Platform:** Illumina HumanHT-12 V4.0 (GPL10558)
**Source:** Grinchuk et al., 2018

| Sample Type | n |
|---|---|
| HCC Tumor | 115 |
| Adjacent Non-tumor | 52 |
| **Total** | **167** |

- Staging: ✅ **BCLC staging available** (in series matrix file)
- Same Illumina platform as GSE36376 — consistent chemistry

---

## Total Replication Coverage
- **3 independent GEO cohorts**: 1,088 additional samples
- **3 staging-capable cohorts** for stage-stratified validation: TCGA-LIHC + GSE14520 + GSE76427
- **2 platforms** beyond TCGA RNA-seq: Affymetrix + Illumina

---

## Analysis Plan Mapping

| Analysis | Primary | Replication |
|---|---|---|
| DGE (tumor vs normal) | TCGA-LIHC | GSE14520, GSE36376, GSE76427 |
| Stage-stratified DGE (early vs late HCC) | TCGA-LIHC | GSE14520 (TNM), GSE76427 (BCLC) |
| DNA Methylation | TCGA-LIHC | N/A (no methylation in GEO sets) |
| Somatic mutation landscape | TCGA-LIHC | N/A |
| Survival analysis (Cox + KM) | TCGA-LIHC | N/A |
| Multi-omics signature | TCGA-LIHC | Validation in GEO if survival data present |

---

## Data Access Notes
- **TCGA-LIHC**: Open access via `TCGAbiolinks` R package
- **GEO datasets**: Open access via `GEOquery` R package
- **Matched normals**: Available on GDC under Sample Type = "Solid Tissue Normal"
- **Germline SNPs**: NOT directly accessible (TCGA germline VCFs require dbGaP)
  → Substituting with in-silico variant annotation from gnomAD, 1000 Genomes, dbSNP using VEP/ANNOVAR

---

## Local file structure (gitignored, not on GitHub)

data/
└── raw/
├── tcga_lihc_rnaseq_se.rds # ~800 MB - SummarizedExperiment, 60660 × 421
├── tcga_lihc_clinical.rds # ~300 KB - clinical data, 377 rows
├── tcga_lihc_meth450_se.rds # PENDING - to download Week 6
└── GDCdata/ # raw downloaded files from GDC

**Confirmed counts (28 Apr 2026 via TCGAbiolinks):**
- Primary Tumor: 371
- Solid Tissue Normal: 50
- Total RNA-seq samples: 421
- Genes (Ensembl IDs): 60,660
- Clinical records: 377

---
## Pending Items
- [ ] Verify GSE14520 supplementary clinical file has all 488 samples annotated
- [ ] Check TCPA portal for phospho-SMAD3 antibody on TCGA-LIHC RPPA panel
- [ ] Confirm AJCC stage breakdown counts for TCGA-LIHC (Stage I/II/III/IV)