# ASD-SZ-WGS-Japan

Code repository for **"Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia"** (Kushima et al., *Nature Communications*, 2026).

## Overview

This repository organises the analysis pipelines used to test two primary questions in a harmonised whole-genome sequencing (WGS) study of 11,335 Japanese participants (597 autism spectrum disorder, 762 schizophrenia, 8,610 population-based controls, 1,366 family members):

1. Whether rare tandem repeat expansion (TRE) burden differs across autism spectrum disorder (ASD), schizophrenia (SZ), and matched controls
2. Whether rare deletions disrupt cell-type-resolved developmental 3D genome (TAD) boundaries differentially across ASD and SZ
3. Whether contextual genetic layers (rare coding/CNV burden, polygenic risk score portability) and a joint-layer integration support a layered cross-disorder architecture

## Repository structure

| Module | Directory | Description |
|--------|-----------|-------------|
| **TRE — Tandem Repeat Expansion** | [`tre/`](tre/) | TRE outlier calling (STRling and ExpansionHunter Denovo) with 5-fold cross-fitted rare-burden regression and cross-caller comparisons (Methods §8). See [`tre/README.md`](tre/README.md) for execution order and module details. |
| **TAD — TAD Boundary Disruption** | [`tad/`](tad/) | Cell-type-resolved TAD boundary disruption analysis with a two-tier endpoint structure: an architecture-level primary endpoint (aggregate Diff_any boundary-bin four-test family ASD/SZ × DEL/DUP plus aggregate ASD-versus-SZ contrast and four-layer aggregate-signal robustness) and a secondary class-level localisation (per-class B′-Firth logistic regression with multinomial heterogeneity testing); shared upstream pipeline (`tad/01_heffel_boundary_master`–`tad/05_wgs_sample_burden`); sensitivity layers (matched-static resampling, Diff_any-versus-Static specificity, exon-exclusion two-test in WGS and MSSNG); external replication and three-cohort inverse-variance-weighted meta-analysis (arrayCGH, MSSNG); and Diff_any-bounded TAD-domain constraint enrichment (Methods §9–11). See [`tad/README.md`](tad/README.md) for the full pipeline order and endpoint mapping. |
| **Joint-layer integration** | [`joint_layer/`](joint_layer/) | Integrative joint-layer logistic regression for ASD and SZ across six progressive risk layers (TAD / TRE / PRS / PTV / CNV / NAHR-CNV) with added-last likelihood-ratio testing, plus a multinomial logistic regression with shared controls for the formal ASD-versus-SZ TAD coefficient contrast (Figure 5; Supplementary Table 11). See [`joint_layer/README.md`](joint_layer/README.md). |

## Documentation

- [`tre/README.md`](tre/README.md) — TRE pipeline (STRling, EHdn, cross-caller); execution order, wrapper→helper correspondence, runtime estimates
- [`tad/README.md`](tad/README.md) — TAD pipeline (17 numbered modules); two-tier endpoint hierarchy, pipeline order, endpoint→figure/table mapping
- [`joint_layer/README.md`](joint_layer/README.md) — Joint-layer integrative logistic regression; layer definitions, output mapping to Figure 5 and Supplementary Table 11, verified numerical results
- [`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md) — Draft text for the manuscript Code availability statement (multiple lengths)
- [`LICENSE`](LICENSE) — MIT License

## Environment

```bash
# Python dependencies
pip install -r requirements.txt

# Or full conda environment (Python + bioinformatics tools)
conda env create -f environment.yml
conda activate tre-burden
```

External tools required on `$PATH`: `samtools >= 1.17`, `bedtools >= 2.31`, `STRling 0.5.2`, `ExpansionHunter Denovo 0.9`, `R >= 4.2` (with `readr`, `dplyr`, `tidyr` for the joint-layer module), plus SLURM for job scheduling.

## Reproducibility note

The wrapper scripts in each module perform **no analytical computation**: they only define execution order, submit helper scripts via `sbatch`, chain jobs with `--dependency=afterok`, and record submitted Job IDs and timestamps in manifest files. The analytical implementation lives in each module's `helpers/` directory (or numbered subdirectories for the TAD module). The joint-layer module is implemented as a single self-contained R script (`joint_layer/asd_sz_layered_logistic_v6.R`).

Reproducing the full analyses additionally requires access to cohort-specific input data and reference resources; data access policies are described in the manuscript's Data availability statement.

## Citation

If you use code from this repository, please cite:

> Kushima et al. *Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia.* Nature Communications (2026).

## License

MIT — see [`LICENSE`](LICENSE).
