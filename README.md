# Multi-Omics Analysis of TGF-β Pathway Dysregulation in HCC

MS Bioinformatics thesis project (COMSATS University Islamabad, 2026–2027).

## Objective
Characterize TGF-β signaling pathway dysregulation in hepatocellular carcinoma across multiple molecular layers — gene expression, DNA methylation, somatic mutations, and regulatory variants — using integrated multi-omics analysis of TCGA-LIHC and GEO cohorts.

## Datasets
See [01_data/inventory.md](01_data/inventory.md) for complete data inventory.
- TCGA-LIHC (n=462): primary cohort
- GSE14520, GSE36376, GSE76427: replication cohorts (n=1,088 combined)

## Pipeline (in development)
1. Snakemake pipeline for preprocessing
2. DESeq2/limma differential expression
3. minfi/ChAMP methylation analysis
4. VEP/ANNOVAR variant annotation
5. clusterProfiler pathway enrichment
6. Cox regression + glmnet survival signature
7. R Shiny dashboard for interactive exploration

## Repository structure
- `01_data/` — data inventory and metadata
- `02_scripts/` — R, Python, and Snakemake code
- `03_results/` — output tables (DEGs, methylation, signatures)
- `04_figures/` — publication-ready figures

## Author
**Tooba Zahra** — MS Bioinformatics, COMSATS University Islamabad
- LinkedIn: https://www.linkedin.com/in/tooba-zahra-ab2015246/
- GitHub: https://github.com/ToobaZahra
- Email: tooba.zahra19@gmail.com

## Supervisor
Dr. Farhan Haq, Tenured Associate Professor, Department of Biosciences

## License
MIT (to be added)