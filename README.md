# str_03242026 — Publication-Ready TRE Burden Analysis Pipeline

## Overview

This directory contains the publication-ready organization of the tandem repeat
expansion (TRE) burden analysis pipeline for the GRIFIN-PD whole-genome
sequencing study.

Two project layers currently coexist:

| Directory | Role |
|-----------|------|
| `str_12282025` | Original working directory containing the helper scripts, intermediate outputs, and full analysis history. This directory is retained unchanged for provenance and reproducibility. |
| `str_03242026` | Reorganized entry-point directory containing clean wrapper scripts, shared symlinked inputs, and this README. This directory is intended for code-availability presentation and reproducible workflow entry points. |

At the current stage, the wrapper scripts in `str_03242026` call helper scripts
that still reside in `str_12282025`. This is intentional. The wrappers improve
workflow readability and execution-order clarity, while the original helper
scripts and outputs remain preserved in the historical analysis directory. A
future step may migrate helper scripts into `str_03242026` for full
self-containment.

## Directory structure

```text
str_03242026/
├── config_v1.sh
├── README.md
├── common/
│   ├── sample_lists/
│   │   ├── ehdn_all_samples.tsv
│   │   └── casecontrol_samples.tsv
│   ├── depth/
│   │   └── depths_all.tsv
│   ├── resources/
│   │   ├── gene_regions_1kb_pad.bed
│   │   └── blacklist/
│   └── logs/
├── ehdn/
│   ├── 01_ehdn_setup_and_sample_lists.sh
│   ├── 02_ehdn_depth.sh
│   ├── 03_ehdn_profile_and_merge.sh
│   ├── 04_ehdn_casecontrol_burden.sh
│   └── logs/
├── strling/
│   ├── 01_strling_build_panel.sh
│   ├── 02_strling_casecontrol_call_and_outliers.sh
│   ├── 03_strling_casecontrol_burden.sh
│   └── logs/
└── crosscaller/
    └── 04_tre_crosscaller_compare.sh
```

### Shared inputs

The `common/` directory contains shared inputs linked from the original
analysis directory:

- `common/sample_lists/` — `ehdn_all_samples.tsv`, `casecontrol_samples.tsv`
- `common/depth/` — `depths_all.tsv`
- `common/resources/` — `gene_regions_1kb_pad.bed` and blacklist/repeat-filtering resources used by EHdn

These symlinks keep the reorganized workflow lightweight while preserving the
original data layout.

## Shared configuration

The helper scripts use two configuration layers depending on the caller.

- **STRling helpers** source `str_12282025/strling/00_strling_genomewide_config_v1.sh`, which defines STRling-specific paths, panel-construction parameters, sample lists, output directories, and SLURM defaults.
- **EHdn helpers** source `str_12282025/config_v1.sh`, which defines project-wide paths, reference genome, EHdn binary, sample metadata, depth files, shared resources, and SLURM defaults.

In the current intermediate layout, the entry-point wrappers reside under
`str_03242026`, whereas the helper scripts and their configuration files still
reside under `str_12282025`.

---

## STRling pipeline — execution order

Run the following wrapper scripts in order. Each wrapper submits downstream
jobs with `--dependency=afterok` so that later steps wait for earlier ones to
finish successfully.

### 1. `strling/01_strling_build_panel.sh`

Build the genome-wide STRling panel and derive the genic 1 kb-padded 3–8 bp
in-bounds panel used in the burden analysis.

### 2. `strling/02_strling_casecontrol_call_and_outliers.sh`

Run STRling genotyping on case-control samples, filter calls to in-bounds
loci, and generate `STRs.tsv` using `strling-outliers.py`.

### 3. `strling/03_strling_casecontrol_burden.sh`

Run the 5-fold cross-fitted rare TRE burden analysis for STRling and then
perform downstream QC and sensitivity analyses.

### 4. `crosscaller/04_tre_crosscaller_compare.sh`

Run the integrated EHdn–STRling ASD-vs-SZ case-case comparison and
heterogeneity analyses.

### STRling wrapper → helper correspondence

| Wrapper | Helpers called |
|---------|---------------|
| `strling/01_strling_build_panel.sh` | `strling/01_strling_extract_array_genomewide_v3.sh` → `strling/03_strling_merge_by_chrom_genomewide_v3.sh` → `strling/04_strling_make_joint_bounds_genomewide_v1.sh` → `strling/04b_make_joint_bounds_genic_len3_8_v2.sh` |
| `strling/02_strling_casecontrol_call_and_outliers.sh` | `strling/06_strling_call_array_genic_v1.sh` → `strling/10_make_calls_genic_inbounds_v1.py` → `strling/11_strling_outliers_casecontrol_inbounds_v1.sh` |
| `strling/03_strling_casecontrol_burden.sh` | `strling/08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py` → `strling/09_strling_qc_sensitivity_rare_inbounds_v4.py` |
| `crosscaller/04_tre_crosscaller_compare.sh` | `ehdn-strling_01022026/20_tre_case_case_comparison_v3.py` |

All helper paths are relative to `str_12282025/`.

---

## EHdn pipeline — execution order

Run the following wrapper scripts in order.

### 1. `ehdn/01_ehdn_setup_and_sample_lists.sh`

Create required directories, verify key resources, and generate
purpose-specific sample lists from the sample metadata.

### 2. `ehdn/02_ehdn_depth.sh`

Compute per-sample average genome depth and aggregate all depth estimates into
`depths_all.tsv`. This depth file is shared across TRE analyses.

### 3. `ehdn/03_ehdn_profile_and_merge.sh`

Run EHdn profiling on all samples, including `family_member` samples, and then
merge the resulting `locus.tsv` files into a genic, depth-normalized matrix.

### 4. `ehdn/04_ehdn_casecontrol_burden.sh`

Run the 5-fold cross-fitted rare TRE burden analysis for EHdn, followed by
downstream burden statistics including logistic regression for carrier status,
Poisson regression with `offset=log(observed_clusters_total)` for event counts,
and Mann–Whitney U tests.

### EHdn wrapper → helper correspondence

| Wrapper | Helpers called |
|---------|---------------|
| `ehdn/01_ehdn_setup_and_sample_lists.sh` | `00_setup_project_v4.sh` → `ehdn/01_prepare_sample_lists_v2.py` |
| `ehdn/02_ehdn_depth.sh` | `strling/01_calc_depth_array_fast_v2.sh` → `strling/03_collect_depths_v1.py` |
| `ehdn/03_ehdn_profile_and_merge.sh` | `ehdn/02_run_ehdn_array_v2.sh` → `ehdn/04_merge_ehdn_novel_norm_v2.py` |
| `ehdn/04_ehdn_casecontrol_burden.sh` | `ehdn/14_outlier_burden_rare_casecontrol_crossfit_v19.py` → `ehdn/17_burden_statistical_test_v20.py` |

All helper paths are relative to `str_12282025/`.

---

## Utilities excluded from the main analysis path

The following scripts are retained for reference but are not part of the
current primary TRE burden workflow:

- `str_12282025/ehdn/15_annotate_genes_v18_gencode.py`
- `str_12282025/ehdn/16_gene_significance_test_v18.py`

These scripts were used in earlier exploratory or legacy analyses and are not
required to reproduce the main burden results reported in the manuscript.

## Result invariance

The wrapper scripts perform **no analytical computation**. They only:

1. determine array sizes from sample lists,
2. submit helper scripts via `sbatch`,
3. chain jobs with `--dependency=afterok`, and
4. record Job IDs and timestamps in manifest files.

**Given identical inputs and identical helper scripts, the analytical results
are unchanged.** This reorganization is intended only to improve
execution-order clarity and code-availability presentation. This applies
equally to both the STRling and EHdn pipelines.

## Logging and manifests

Each wrapper writes:

- a timestamped manifest file under its own `logs/` directory,
- submitted SLURM Job IDs, and
- wrapper-level stdout/stderr logs where applicable.

The underlying helper scripts continue to write their original logs in the
historical analysis directory unless and until the helper layer is migrated.

## Approximate runtime

### STRling (9,969 case-control samples, ncbn-cpu partition)

| Step | Approximate wall-clock time | Notes |
|------|-----------------------------|-------|
| `01` / extract | ~1 h | Array job, 40 parallel |
| `01` / merge | ~8 h | Longest chromosome typically dominates |
| `01` / joint-bounds | ~50 min | Single job |
| `01` / genic filter | ~2 s | Single job |
| `02` / call | ~2–3 min per sample | Array job, 50 parallel |
| `02` / in-bounds filter | ~8 s | Single job |
| `02` / outlier detection | ~6 h | Single job, 256 GB RAM |
| `03` / burden | ~2 min | Single job |
| `03` / QC | ~1 s | Single job |
| `04` / case-case comparison | ~1 s | Single job |

### EHdn (11,386 samples, ncbn-cpu partition)

| Step | Approximate wall-clock time | Notes |
|------|-----------------------------|-------|
| `01` / setup | ~1 s | Single job |
| `01` / sample lists | ~30 s | Single job |
| `02` / depth array | ~1 h | Array job, 100 parallel |
| `02` / collect depths | ~1 min | Single job |
| `03` / EHdn profile | ~24–48 h | Array job, 40 parallel, 8 CPUs/task |
| `03` / merge | ~30 min | Single job, 64 GB RAM |
| `04` / burden | ~2–4 h | Single job, 128 GB RAM |
| `04` / burden statistics | ~30 min | Single job |

## Recommended usage

For publication-oriented reproduction, use the wrapper scripts in
`str_03242026` rather than calling helper scripts manually. This provides a
cleaner execution history, explicit job dependency handling, and clearer
code-availability presentation while preserving the original analytical
implementation.
