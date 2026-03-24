# TRE Burden Analysis Pipeline — GRIFIN-PD WGS

## Overview

This repository contains the publication-ready workflow for the tandem repeat
expansion (TRE) burden analysis in the GRIFIN-PD whole-genome sequencing study.

The repository is organized into two layers:

| Layer | Contents |
|-------|----------|
| Entry-point wrappers (`strling/`, `ehdn/`, `crosscaller/`) | Clean wrapper scripts that define execution order, SLURM job dependencies, and logging. These are the recommended starting point for understanding and reproducing the workflow. |
| Helper scripts (`helpers/`) | The underlying analytical scripts called by the wrappers. These perform the actual computation (genotyping, merging, burden testing, etc.). |

The wrapper scripts perform **no analytical computation**. They only determine
array sizes from sample lists, submit helper scripts via `sbatch`, chain jobs
with `--dependency=afterok`, and record Job IDs and timestamps in manifest
files.

## Directory structure

```text
.
├── config_v1.sh                          # Central configuration
├── README.md
├── CODE_AVAILABILITY.md
├── .gitignore
├── strling/                              # STRling entry-point wrappers
│   ├── 01_strling_build_panel.sh
│   ├── 02_strling_casecontrol_call_and_outliers.sh
│   └── 03_strling_casecontrol_burden.sh
├── ehdn/                                 # EHdn entry-point wrappers
│   ├── 01_ehdn_setup_and_sample_lists.sh
│   ├── 02_ehdn_depth.sh
│   ├── 03_ehdn_profile_and_merge.sh
│   └── 04_ehdn_casecontrol_burden.sh
├── crosscaller/                          # Cross-caller entry-point wrapper
│   └── 04_tre_crosscaller_compare.sh
└── helpers/                              # Underlying analytical scripts
    ├── 00_setup_project_v4.sh
    ├── strling/
    │   ├── 00_strling_genomewide_config_v1.sh
    │   ├── 01_strling_extract_array_genomewide_v3.sh
    │   ├── 01_calc_depth_array_fast_v2.sh
    │   ├── 03_strling_merge_by_chrom_genomewide_v3.sh
    │   ├── 03_collect_depths_v1.py
    │   ├── 04_strling_make_joint_bounds_genomewide_v1.sh
    │   ├── 04b_make_joint_bounds_genic_len3_8_v2.sh
    │   ├── 06_strling_call_array_genic_v1.sh
    │   ├── 08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py
    │   ├── 09_strling_qc_sensitivity_rare_inbounds_v4.py
    │   ├── 10_make_calls_genic_inbounds_v1.py
    │   └── 11_strling_outliers_casecontrol_inbounds_v1.sh
    ├── ehdn/
    │   ├── 01_prepare_sample_lists_v2.py
    │   ├── 02_run_ehdn_array_v2.sh
    │   ├── 04_merge_ehdn_novel_norm_v2.py
    │   ├── 14_outlier_burden_rare_casecontrol_crossfit_v19.py
    │   └── 17_burden_statistical_test_v20.py
    └── crosscaller/
        └── 20_tre_case_case_comparison_v3.py
```

## Path configuration

The helper scripts contain hardcoded paths from the original analysis
environment (NIG supercomputer). These paths are marked with `# CONFIGURE`
comments. To run the pipeline in a different environment, update the following:

- **`config_v1.sh`**: Central configuration file. Paths can be overridden via
  environment variables using the `${VAR:-default}` pattern.
- **`helpers/strling/00_strling_genomewide_config_v1.sh`**: STRling-specific
  configuration (panel paths, output directories, SLURM defaults).
- **Individual helper scripts**: Some scripts define `PROJECT_ROOT` and
  `OUT_ROOT` internally; update these or set them as environment variables.

Key paths that need adjustment:

| Variable | Description | Original value |
|----------|-------------|----------------|
| `PROJECT_ROOT` | Analysis root directory | `/lustre12/home/kushima-pg/str_12282025` |
| `REFERENCE_FASTA` | GRCh38 reference genome | `/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta` |
| `SAMPLE_INFO` | Sample metadata file | `/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt` |
| `PCA_EIGENVEC` | PCA eigenvector file | `/lustre12/home/kushima-pg/PRS/.../pca.eigenvec` |
| `CRAM_BASE_DIR1` | CRAM directory (GRIFIN) | `/lustre12/home/grifinpd-pg/analysis/parabricks` |
| `CRAM_BASE_DIR2` | CRAM directory (NCBN) | `/lustre12/home/ncbn-share-pg/control_genome/pb3.1.0/results` |

Sample list files (`ehdn_all_samples.tsv`, `casecontrol_samples.tsv`) and
depth files (`depths_all.tsv`) are excluded from this repository because they
contain sample-level identifiers.

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
| `strling/01_strling_build_panel.sh` | `helpers/strling/01_strling_extract_array_genomewide_v3.sh` → `helpers/strling/03_strling_merge_by_chrom_genomewide_v3.sh` → `helpers/strling/04_strling_make_joint_bounds_genomewide_v1.sh` → `helpers/strling/04b_make_joint_bounds_genic_len3_8_v2.sh` |
| `strling/02_strling_casecontrol_call_and_outliers.sh` | `helpers/strling/06_strling_call_array_genic_v1.sh` → `helpers/strling/10_make_calls_genic_inbounds_v1.py` → `helpers/strling/11_strling_outliers_casecontrol_inbounds_v1.sh` |
| `strling/03_strling_casecontrol_burden.sh` | `helpers/strling/08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py` → `helpers/strling/09_strling_qc_sensitivity_rare_inbounds_v4.py` |
| `crosscaller/04_tre_crosscaller_compare.sh` | `helpers/crosscaller/20_tre_case_case_comparison_v3.py` |

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
| `ehdn/01_ehdn_setup_and_sample_lists.sh` | `helpers/00_setup_project_v4.sh` → `helpers/ehdn/01_prepare_sample_lists_v2.py` |
| `ehdn/02_ehdn_depth.sh` | `helpers/strling/01_calc_depth_array_fast_v2.sh` → `helpers/strling/03_collect_depths_v1.py` |
| `ehdn/03_ehdn_profile_and_merge.sh` | `helpers/ehdn/02_run_ehdn_array_v2.sh` → `helpers/ehdn/04_merge_ehdn_novel_norm_v2.py` |
| `ehdn/04_ehdn_casecontrol_burden.sh` | `helpers/ehdn/14_outlier_burden_rare_casecontrol_crossfit_v19.py` → `helpers/ehdn/17_burden_statistical_test_v20.py` |

---

## Utilities excluded from the main analysis path

The following scripts were used in earlier exploratory analyses and are not
part of the primary TRE burden workflow reported in the manuscript:

- `ehdn/15_annotate_genes_v18_gencode.py`
- `ehdn/16_gene_significance_test_v18.py`

## Result invariance

The wrapper scripts perform **no analytical computation**. Given identical
inputs and identical helper scripts, the analytical results are unchanged.
This reorganization is intended only to improve execution-order clarity and
code-availability presentation. This applies equally to both the STRling and
EHdn pipelines.

## Logging and manifests

Each wrapper writes a timestamped manifest file under its own `logs/`
directory, including submitted SLURM Job IDs and wrapper-level timing.

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

For publication-oriented reproduction, use the wrapper scripts in the
top-level `strling/`, `ehdn/`, and `crosscaller/` directories rather than
calling helper scripts directly. This provides a cleaner execution history,
explicit job dependency handling, and clearer workflow presentation while
preserving the original analytical implementation in `helpers/`.
