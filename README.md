# Sorghum Pareto–Phenolics Analysis Pipeline

Reproducible analysis pipeline for the manuscript:

**“Pareto geometry and entropy–yield scaling map trade-offs in sorghum phenolic mutants.”**

This repository rebuilds all derived results from publicly available input tables (Lee et al., 2023) and writes
the combined workbook **`Supplementary_Data_S1_master.xlsx`** plus CSV exports under `./outputs/`.

## What this repo does

Main pipeline (`sorghum_analysis_pipeline.py`) performs:

- Phenotypic embedding (PCA; optional UMAP if installed)
- Multi-objective Pareto optimization (efficient frontier + knee point)
- Differential phenolic co-variation (correlation difference + Fisher’s z-test)
- Compositional entropy (Shannon H) and derived phenotype summaries
- Correlation-graph “control shift” (eigenvector centrality)
- Allele-frequency shift summaries on GWAS-hit SNPs (ΔAF)
- Genetic convergence / equifinality summaries on key SNPs

## Data source

The raw phenotypic and genotypic tables were originally published by:

Lee, Y.-J. et al. (2023). *Genome-wide association study (GWAS) of the agronomic traits and phenolic content in sorghum
(Sorghum bicolor L.) genotypes.* Agronomy 13, 1449. https://doi.org/10.3390/agronomy13061449

This repository does not claim ownership of the original datasets. Please cite the paper above when using the input tables.

## Inputs

Place the following Excel files either in the repository root or in `./data/`:

- `Table S2.xlsx` — phenotypes (includes the 6 phenolic traits; row order defines sample_id 1–96)
- `Table S3.xlsx` — SNP summary metrics per line (includes total SNP count and heterozygous SNP count)
- `Table S4.xlsx` — genotype matrix (SNPs × samples; script canonicalizes sample columns like `23_1 → 23`)
- `Table S5.xlsx` — GWAS hits (part 1)
- `Table S6.xlsx` — GWAS hits (part 2)

## Quick start (recommended: conda)

### 1) Create environment
```bash
conda env create -f environment.yml
conda activate gwas_env
```

### 2) Run the main pipeline
```bash
python sorghum_analysis_pipeline.py
```

Outputs are written to:
- `./outputs/tables/` (CSV tables + a combined Excel workbook)

## Key outputs

- `outputs/tables/Supplementary_Data_S1_master.xlsx`
  - `phenotypes_master`, `pareto_front`
  - `diff_correlation`, `diff_corr_pvalues`
  - `entropy_data`, `control_shift`
  - `hit_snp_pairs`, `hit_snp_altfreq_delta`, `key_snps_for_convergence`
  - `genetic_convergence`, `convergence_stats`
  - `sampleid_to_line`

- `outputs/tables/run_summary.json` (run metadata: Pareto count, knee line, thresholds, seed)

Reproducibility notes:
- Random seed is fixed (42) for stable manifold projection.
- ΔAF uses: alt allele != ref (robust).
- Genetic convergence PCA encodes the second allele in the `Allele` field to match the manuscript behavior.

## Robustness / sensitivity suite (add-on)

Utilities supporting robustness checks and objective/proxy sensitivity analyses:

- `tools/check_s1_reproducibility.py`  
  Compare two `Supplementary_Data_S1_master.xlsx` files (or two `outputs/tables/` directories).

- `tools/sorghum_robustness_suite.py`  
  Writes additional tables under `outputs/tables/reviewer_response/`:
  - entropy–yield correlations (overall + by type)
  - Pareto sensitivity to SNP proxy choice (Jaccard overlap)
  - jackknife stability of Pareto membership and knee selection
  - top-k vs Pareto overlap
  - objective-switch top-5 summaries

- `tools/run_reproduce_and_robustness.py`  
  One-command workflow: rerun pipeline → (optional) compare S1 → run robustness suite.

Example:
```bash
python tools/run_reproduce_and_robustness.py \
  --repo . \
  --data /path/to/folder_with_Table_S2_to_S6 \
  --baseline /path/to/baseline_outputs/tables
```

## License / citation

Please cite the associated manuscript and Lee et al. (2023) when using these scripts.
