# Multi-Omics Analysis of TGF-β Pathway Dysregulation in HCC

MS Bioinformatics thesis project (COMSATS University Islamabad, 2026–2027).

> 🚧 **Status: Active development** — Results will be made public upon thesis submission (January 2027)

---

## Objective
Characterize TGF-β signaling pathway dysregulation in hepatocellular carcinoma across multiple molecular layers — gene expression, DNA methylation, somatic mutations, and regulatory variants — using integrated multi-omics analysis of TCGA-LIHC and GEO cohorts.

---

## Datasets
See [01_data/inventory.md](01_data/inventory.md) for complete data inventory.

| Cohort | Role | n samples |
|--------|------|-----------|
| TCGA-LIHC | Primary analysis | 462 |
| GSE14520 | DGE replication | 488 |
| GSE36376 | DGE replication | 433 |
| GSE76427 | DGE + stage replication | 167 |

---

## Analysis Pipeline

| Step | Method | Status |
|------|--------|--------|
| Data download | TCGAbiolinks, GEOquery | ✅ Done |
| DESeq2 DGE (tumor vs normal) | DESeq2 | ✅ Done |
| Stage-stratified DGE | DESeq2 | ✅ Done |
| GEO replication (3 cohorts) | limma | ✅ Done |
| DNA methylation 450K | limma | ✅ Done |
| Survival analysis | Cox regression + KM | ✅ Done |
| Pathway enrichment | clusterProfiler | ✅ Done |
| Multi-omics integration | Custom R pipeline | ✅ Done |
| Variant annotation | VEP/ANNOVAR | ⏭ Next semester |
| R Shiny dashboard | Shiny | ⏭ Next semester |

---

## Key Findings (preliminary)
- **16/26** TGF-β nexus genes replicated consistently across all 4 cohorts
- **HDAC11** identified as top multi-omics hit: upregulated in tumor, hypomethylated, poor prognosis marker (Cox HR=1.29, p=0.009), replicated in independent cohorts
- Global hypomethylation in HCC tumor (88,102 hypomethylated vs 22,152 hypermethylated probes)
- GSEA confirms proliferation activation and metabolic suppression in HCC

---

## Repository Structure
- `01_data/` — data inventory and metadata
- `02_scripts/` — R analysis scripts (20 scripts, fully documented)
- `03_results/` — output tables (hidden until thesis submission)
- `04_figures/` — publication-ready figures (hidden until thesis submission)

---

## Author
**Tooba Zahra** — MS Bioinformatics, COMSATS University Islamabad
- LinkedIn: https://www.linkedin.com/in/tooba-zahra-ab2015246/
- GitHub: https://github.com/ToobaZahra
- Email: tooba.zahra19@gmail.com

## Supervisor
Dr. Farhan Haq, Tenured Associate Professor, Department of Biosciences

## License
© 2026 Tooba Zahra. Code is available for academic use with attribution (CC BY-NC 4.0).
Results and figures will be released upon thesis submission (January 2027).