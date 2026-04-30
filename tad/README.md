# TAD — TAD Boundary Disruption Pipeline

Code for the TAD boundary disruption analysis used in **"Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia"** (Kushima et al., *Nature Communications*, 2026), corresponding to Methods §9–11 and Figures 3–4 (and related Supplementary Tables / Figures).

This module is a sibling of [`tre/`](../tre/). For the cohort, design, and overall study scope, see the [repository root README](../README.md).

## Overview

The TAD pipeline tests whether rare deletions disrupt cell-type-resolved developmental 3D genome (TAD) boundaries differentially across ASD, SZ, and matched controls. The pipeline produces:

- **Discovery (WGS, Japanese cohort)** — B′ logistic regression at the L2 cell-type level (10 classes; HPC and PFC, Excitatory/Inhibitory/Astro), with carrier-vs-noncarrier sensitivity, an exon-exclusion sensitivity (separate / joint / SV-level), a matched-static resampling, and a Diff_any-vs-Static specificity test.
- **External replication** — arrayCGH (Japanese ASD case-control) and MSSNG (Canadian ASD families) under harmonised covariates and within-cohort BH-FDR.
- **Three-cohort integration** — Inverse-variance-weighted meta-analysis (WGS + arrayCGH + MSSNG).
- **Constraint enrichment** — Diff_any-bounded TAD domain constraint analysis (gnomAD pLI / LOEUF / mis-Z) versus Static-bounded controls.

## Pipeline Modules

Modules are numbered to reflect dependency order (downstream depends on upstream). Each module is a self-contained directory with the analytical script (`*.py` or `*.R`), a SLURM wrapper (`*.sbatch`), and—where relevant—pattern-specific sensitivity wrappers.

| # | Module | Latest version | Paper reference |
|---|--------|----------------|-----------------|
| 01 | `01_heffel_boundary_master/` — Heffel TAD boundary master construction | v9 | Methods §9.1; Supp Table 3 |
| 02 | `02_bin_l2_annotation/` — 25-kb bin × L2 cell-type annotation | **v3 (Q1)** | Methods §9.2 |
| 03 | `03_wgs_sv_events/` — Rare WGS SV (DEL/DUP) extraction | v9 | Methods §9.3 |
| 04 | `04_wgs_sv_boundary_overlap/` — SV × TAD boundary overlap (Pattern A primary + Pattern B/C sensitivity) | v10 / v11 | Methods §9.4 |
| 05 | `05_wgs_sample_burden/` — Per-sample × L2 × SV-type burden table | **v4 (Q1, Pattern A) / v5 (Pattern B/C)** | Methods §9.5 |
| 06 | `06_wgs_primary_L2/` — B′ logistic regression on 10 L2 classes (Pattern A primary + Pattern B/C sensitivity) and MNLogit heterogeneity test | **v5.R (Q1, Pattern A) / v6.R (Pattern B/C)** | Figure 3a; Supp Table 5 |
| 07 | `07_wgs_matched_static/` — Diff_any-vs-Static matched resampling specificity test | v7 / v8 sbatch | Figure 3b; Supp Table 6 |
| 08 | `08_arraycgh_sample_burden/` — arrayCGH replication cohort burden | v22 | Methods §9.6; Supp Table 7 |
| 09 | `09_mssng_sample_burden/` — MSSNG replication cohort burden + meta-analysis | v19 + meta v12 | Methods §9.6 / §10; Supp Table 8 |
| 10 | `10_replication_2way_meta/` — Three-cohort inverse-variance-weighted meta-analysis | v5 | Figure 3c; Supp Table 9 |
| 11 | `11_diff_all_vs_static/` — Diff_any-vs-Static specificity test (logistic regression) | v5 | Figure 3d; Supp Table 10 |
| 12 | `12_exon_exclusion_wgs/` — WGS exon-exclusion sensitivity (A1/A2/A3/A3-strict) | **v4 (Q2)** | Figure 4a; Supp Table 11 |
| 13 | `13_exon_exclusion_mssng/` — MSSNG exon-exclusion sensitivity | **v4 (Q2)** | Figure 4b; Supp Table 12 |
| 14 | `14_constraint_enrichment_v1/` — Diff_any-bounded TAD domain constraint enrichment (Diff_any vs Static pipeline) | v9 / static_v1 | Figure 4c–d; Supp Table 13 |
| 15 | `15_l2_jaccard_v1/` — L2 cell-type Jaccard similarity across boundaries | v2 | Supp Figure / Supp Methods |

**Q1 / Q2 / Q3 update annotations** (April 2026, repository housekeeping; see "Reproducibility note" below):
- **Q1**: Pattern A primary scripts (Modules 02 / 05 / 06) had paper-irrelevant S2–S5 sensitivity grouping removed; numerical outputs are bit-identical to the paper-canonical baseline.
- **Q2**: Modules 12 + 13 share a common helper module (`common/exon_helpers.py`); ~310 lines of duplicate code were factored out without altering numerical outputs.
- **Q3**: All path constants are routed through `common/paths_v1.py` with environment-variable overrides (NIG NCBN defaults are preserved); see `config_v1.sh`.

## Pattern A primary vs Pattern B / C sensitivity

The TAD pipeline runs three SV filtering patterns in parallel for sensitivity:

| Pattern | NAHR-mediated SV | Use |
|---------|-------------------|-----|
| **Pattern A (primary)** | excluded | Main results in Figure 3 / Figure 4 |
| Pattern B (sensitivity) | included | Robustness check (Supp Tables) |
| Pattern C (sensitivity) | included with size cutoff | Robustness check (Supp Tables) |

The wrapper scripts `*_v11_patB.sbatch` (Module 04), `*_v5_patB.sbatch` (Module 05), `*_v6_patB.sbatch` (Module 06) execute Pattern B; `*_patC.sbatch` execute Pattern C. Pattern A primary uses the unsuffixed wrappers. The same analytical scripts (e.g. `34_intersect_sv_with_heffel_master_v11.py`, `38_compute_sample_burden_L2_and_specificity_v5.py`, `39_fit_B_prime_L2_and_specificity_v6.R`) drive Pattern B and C — only the SV filter input differs.

## Common helpers

| File | Purpose |
|------|---------|
| `common/paths_v1.py` | Single source of truth for all input/output paths. Resolves `PIPELINE_ROOT` from environment variable (defaults to NIG NCBN paths). |
| `common/paths_v1_manifest.json` | Path manifest in JSON for R-side import (R does not source `paths_v1.py` directly). |
| `common/exon_helpers.py` | Shared helpers for Modules 12 + 13 (GENCODE exon merge, bin-exon overlap, L2 annotation loader). |
| `common/naming_v1.py` | L2 class / SV-type / comparison label canonicalisation. |
| `config_v1.sh` | Environment-variable exporter; source this before running `sbatch` if you need to override defaults. |

## Reproducibility verification

`99_verify_vs_draft/verify_paper_numbers_v3.py` empirically compares the Module 06 v5 output against the primary numbers reported in the manuscript:

- HPC_Exc-DG (ASD vs HC, n_boundary, DEL): OR = 2.79, 95% CI 1.59–4.88, BH-FDR = 1.7 × 10⁻³ (strongest L2 class)
- 9 / 10 L2 classes BH-FDR < 0.05 (PFC_Astro the only non-significant class)
- OR range across the 9 significant classes: 1.77 – 2.79
- SZ vs HC and DUP: no significant positive class

Run from the repo root after sourcing `config_v1.sh`:

```bash
python3 99_verify_vs_draft/verify_paper_numbers_v3.py
# Expected output: PASS: 10, FAIL: 0
```

## Environment

The analytical scripts assume the following on `$PATH`:

- `python3 >= 3.10` with `pandas`, `numpy`, `scipy`, `statsmodels` (with GEE), `tqdm`
- `Rscript >= 4.2` with `data.table`, `dplyr`, `tidyr`, `purrr`
- `bedtools >= 2.31`
- SLURM (`sbatch`) for job scheduling

Reference resources used (paths defined in `common/paths_v1.py`):

- GENCODE v46 protein-coding exon GTF (`/lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz`)
- Heffel et al. (2024) HPC + PFC TAD boundary master (preprocessed to `01_heffel_boundary_master/output_v9/heffel_boundary_master_v9.tsv.gz`)

## Reproducibility note

Running this pipeline end-to-end requires access to:

1. The NIG NCBN cohort WGS cram files (controlled access; see manuscript's Data availability statement).
2. The arrayCGH replication cohort SV table (Japanese ASD case-control; see manuscript).
3. The MSSNG replication cohort SV calls (DECIPHER / DCC controlled access; see manuscript).
4. The Heffel et al. (2024) TAD boundary HPC + PFC master tables (downloadable from the original publication's data release).

The default paths in `common/paths_v1.py` reflect the NIG NCBN compute environment used during manuscript preparation. Adapt these (or set the corresponding environment variables in `config_v1.sh`) for use on a different system.

## Versioning

This module reflects the state of `tad04292026/` on the NIG NCBN cluster as of April 30, 2026. Earlier internal iterations (v1 – v8 etc.) are not committed; only the production (latest) version of each script is included.

## Citation

If you use code from this module, please cite:

> Kushima et al. *Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia.* Nature Communications (2026).

## License

MIT — see [`../LICENSE`](../LICENSE) at the repository root.
