# Sorghum Pareto–Phenolics Analysis Pipeline

This repository provides a reproducible analysis pipeline for the manuscript:

**“An efficient frontier links phenolic specialization to high-yield states in mutagenized sorghum.”**

The pipeline performs:
- Phenotypic manifold learning (PCA + optional UMAP)
- Multi-objective Pareto optimization (efficient frontier and knee point)
- Differential metabolic network rewiring (correlation difference + Fisher’s Z-test)
- Metabolic specialization (Shannon entropy of phenolic composition)
- Metabolic control shift (eigenvector centrality)
- Genetic convergence analysis on elite-associated loci (PCA on key GWAS loci)

## Data provenance (source of input files)

The **raw phenotypic and genotypic datasets** analyzed here were originally published by:

Lee, Y.-J. et al. (2023). *Genome-wide association study (GWAS) of the agronomic traits and phenolic content in sorghum (Sorghum bicolor L.) genotypes.* **Agronomy** 13, 1449. https://doi.org/10.3390/agronomy13061449

This repository does **not** claim ownership of the original datasets. Please cite the paper above when using the input tables.

## Inputs

Place the following files either in the repository root or in `./data/`:

- `Table S2.xlsx` — Phenotypes (includes the 6 phenolic traits used for the manifold)
- `Table S3.xlsx` — SNP burden per line
- `Table S4.xlsx` — Genotypes (SNP matrix; includes sample columns 1–96; the script canonicalizes columns like `23_1` → `23`)
- `Table S5.xlsx` — GWAS hits (part 1)
- `Table S6.xlsx` — GWAS hits (part 2)

> Note: Users may obtain these as supplementary tables from Lee et al. (2023) and/or associated public archives referenced in that publication.

## Quick start

### 1) Install dependencies
```bash
pip install -r requirements.txt


### 2) Run the pipeline
python sorghum_analysis_pipeline.py

utputs will be created under:

./outputs/tables/ (CSV tables + a combined Excel file)

Outputs

Key outputs include:

outputs/tables/Supplementary_Data_S1_master.xlsx
A combined workbook containing:

phenotypes_master, pareto_front

diff_correlation, diff_corr_pvalues

entropy_data, control_shift

hit_snp_pairs, hit_snp_altfreq_delta, key_snps_for_convergence

genetic_convergence, convergence_stats

sampleid_to_line

CSV exports of the same tables for convenient downstream use.

Reproducibility notes

Random seed is fixed (42) to ensure stable manifold projection.

Delta allele frequency (ΔAF) uses alt allele ≠ ref (robust).

Genetic convergence PCA encodes the second allele (index 1) in the Allele field to match the manuscript’s reported behavior.

Data availability statement (for the manuscript)

The raw phenotypic and genotypic datasets analyzed in this study were originally published by Lee et al. (2023) and are publicly available.
All derived data generated in this study, including Pareto optimality labels, metabolic entropy metrics, and network analysis results, are provided in Supplementary Data 1 accompanying the manuscript.
All custom Python scripts used to reproduce the analyses are available in this repository.
