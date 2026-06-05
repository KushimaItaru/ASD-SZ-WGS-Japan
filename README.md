# ASD-SZ-WGS-Japan

[![DOI](https://zenodo.org/badge/1164528445.svg)](https://doi.org/10.5281/zenodo.20556335)

Code repository for **"Shared rare tandem-repeat expansion burden and autism-enriched developmental TAD-boundary-annotated deletion burden across autism and schizophrenia"** (Kushima et al., 2026; manuscript under review).

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

External tools required on `$PATH`: `samtools >= 1.17`, `bedtools >= 2.31`, `STRling 0.5.2`, `ExpansionHunter Denovo 0.9`, `R >= 4.2` (`readr`, `tibble`, `tidyr`, `dplyr`, `stringr`, `purrr`, `jsonlite`, `ggplot2`, `data.table`, `scales`, `pROC`), plus SLURM for job scheduling.

## Reproducibility note

The wrapper scripts in each module perform **no analytical computation**: they only define execution order, submit helper scripts via `sbatch`, chain jobs with `--dependency=afterok`, and record submitted Job IDs and timestamps in manifest files. The analytical implementation lives in each module's `helpers/` directory (or numbered subdirectories for the TAD module). The joint-layer module is implemented as a single self-contained R script (`joint_layer/asd_sz_layered_logistic_v6.R`).

Reproducing the full analyses additionally requires access to cohort-specific input data and reference resources; data access policies are described in the manuscript's Data availability statement.

## Citation

If you use code from this repository, please cite:

> Kushima et al. *Shared rare tandem-repeat expansion burden and autism-enriched developmental TAD-boundary-annotated deletion burden across autism and schizophrenia.* 2026 (manuscript under review).

## License

MIT — see [`LICENSE`](LICENSE).
<!-- Paste this section into the repository README.md (e.g., after "Environment"). -->

## Demo

A small, fully self-contained demo is provided in [`demo/`](demo/) so that editors and
reviewers can run the primary burden model **without access to the controlled-access
cohort data**.

### What it does
`demo/run_demo_burden.py` fits the B′-aligned logistic burden regression used for the
primary case–control analyses

```
case ~ burden_count + sex + PC1–PC10 + log1p(total_bases) + log1p(total_gene)
```

on a small **synthetic** per-individual table (`demo/demo_burden_input.tsv`; 600
simulated individuals, no real participant data). A modest case-enriched burden is
planted so the model returns a clear positive estimate.

### System requirements
- Operating system: Linux / macOS / Windows (any OS with Python).
- Python ≥ 3.10 with `pandas`, `numpy`, `statsmodels` (see `requirements.txt` /
  `environment.yml`). Tested with Python 3.13, statsmodels 0.14.4.
- No non-standard hardware. Typical install time (dependencies): ~5–15 min;
  the demo itself needs no compilation.

### Run
```bash
cd demo
python make_demo_data.py      # (optional) regenerate the synthetic table (seed = 42)
python run_demo_burden.py
```

### Expected output (seed-fixed; runs in < 5 s on a normal desktop)
```
N (controls / cases) : 400 / 200
  Odds ratio (OR)    : 1.555
  95% CI             : 1.322 - 1.829
  P value            : 9.548e-08
Total run time       : ~0.01 s
```
(Exact values are reproducible with the provided seed; minor differences may occur
across library versions.)

### Notes
- The demo data are **simulated**: they illustrate the analytical method only and do
  not reproduce any manuscript result.
- The full analyses require controlled-access cohort data (see the manuscript Data
  availability statement). Numerical reproducibility of the **reported** numbers is
  provided separately by `tad/99_verify_vs_draft/`, which re-derives the primary
  burden values from intermediate, non-identifying summary inputs.
