# Cell-type-resolved TAD boundary disruption analysis

This directory contains the analysis pipeline for the cell-type-resolved 3D-genome TAD boundary disruption analysis described in *Kushima et al., Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia*.

The pipeline tests whether rare deletions disrupting cell-type-differential TAD-boundary bins are enriched in ASD or SZ relative to controls, using the Heffel et al. developmental 3D genome atlas (snm3C-seq Level-2 [L2] differential boundaries across 11 prefrontal/hippocampal cell-type classes).

## Endpoint hierarchy (matches paper Methods §12 and Results §22–28)

The TAD analysis uses a **two-tier endpoint structure** with the aggregate Diff_any union as the architecture-level **primary endpoint** and per-class L2 results as **secondary localisation**.

### Primary endpoint (architecture-level)

The primary endpoint is the aggregate Diff_any boundary-bin burden, defined as the union of all 25-kb bins flagged as differential in at least one of the ten analysed L2 cell-type classes (4,980 unique bins). Four pre-specified primary tests were performed: ASD-versus-Controls and SZ-versus-Controls deletion (DEL) and duplication (DUP) burden under a B′-aligned logistic regression with Sex, PC1–PC10, log1p(total_bases), and log1p(total_gene_overlap) covariates. Hierarchical BH-FDR was applied within the four-test family (q = 1.2 × 10⁻⁴ for the headline ASD-DEL signal). Aggregate-level ASD-versus-SZ contrast (1 test) and 4-layer aggregate-signal robustness (standard LOCO + strict LOCO + top-N high-burden sample exclusion + top-N high-burden event exclusion) are also part of the primary inference layer.

→ Implemented in **`16_aggregate_primary_sensitivity_v1/`** (computes Fig. 3a Block 1 and Block 3, paper §22, §26, §27).

### Secondary endpoint (class-level localisation)

The secondary endpoint reports per-class L2 results (10 cell-type classes after exclusion of HPC Astro for boundary-bin count < 500) using the same B′ logistic regression model. Class-level results are interpreted as localisation of the primary aggregate signal rather than as independent inference. BH-FDR is applied within the 10-class family. Per-class ASD-versus-SZ heterogeneity is assessed by multinomial logistic regression with shared controls.

→ Implemented in **`06_wgs_primary_L2/`** (computes Fig. 3b, paper §22 supporting; multinomial heterogeneity in `40_fit_MNLogit_heterogeneity_v3.py`, paper §27).

### Sensitivity layers (full list)

- **Matched-static specificity** (`07_wgs_matched_static/tad_dynamic_boundary_specificity_test_v7.py`): per-class observed-OR vs property-matched static-boundary null distribution (10,000 resamples; Fig. 3c, paper §24).
- **Diff_any vs Static per-bin specificity** (`11_diff_all_vs_static/42_compute_diff_any_vs_static_v5.py`): per-bin OR contrast within the same B′ covariate framework (Fig. 3d, paper §24).
- **External replication / IVW meta-analysis** (`08_arraycgh_sample_burden`, `09_mssng_sample_burden`, `10_replication_2way_meta`): Japanese arrayCGH cohort within-array IVW, MSSNG GEE family-clustered, two-cohort IVW meta-analysis (Fig. 3a Block 2 + Fig. 3b external panel + Fig. 3f primary three-cohort IVW for constraint enrichment).
- **Constraint-domain enrichment** (`14_constraint_enrichment_v1/`): three-cohort IVW combining WGS Discovery, arrayCGH, and MSSNG over Diff_any-bounded TAD-scale regulatory domains (Fig. 3f, paper §29).
- **Exon-exclusion** (`12_exon_exclusion_wgs/`, `13_exon_exclusion_mssng/`): A1 primary (exon-free boundary bins) and A2 joint co-exposure analyses in discovery WGS and MSSNG (paper §25). arrayCGH was not used as a formal exon-exclusion replication platform because clinical aCGH probe designs preferentially cover coding exons and have limited sensitivity for exon-poor and intergenic CNVs.
- **arrayCGH cross-disorder SZ vs Controls** (implemented within `08_arraycgh_sample_burden/tad_replication_arraycgh_v22.py`, Phases 6B–6C): exploratory SZ-versus-Control deletion burden + ASD-versus-SZ heterogeneity in the arrayCGH cohort (paper §27 supporting).
- **Reproducibility verification** (`99_verify_vs_draft/`): scripts that re-derive the primary L2 burden numbers reported in the manuscript text against draft-stage reference outputs.

## Pipeline order (numbered steps)

| # | Directory | Purpose | Endpoint role | Paper reference |
|---|---|---|---|---|
| 01 | `01_heffel_boundary_master/` | Build Heffel L2 boundary master (11 classes; differential vs static) | shared upstream | Methods §12.1 |
| 02 | `02_bin_l2_annotation/` | Annotate 25-kb bins with per-class L2 differential / static / Diff_any union | shared upstream | Methods §12.1, §12.6 |
| 03 | `03_wgs_sv_events/` | Extract rare DEL/DUP events from WGS Discovery (Pattern A/B/C SV-filter configurations) | shared upstream | Methods §12.2 |
| 04 | `04_wgs_sv_boundary_overlap/` | Intersect rare SV events with L2 boundary bins (BEDTools) | shared upstream | Methods §12.3 |
| 05 | `05_wgs_sample_burden/` | Per-individual boundary-bin burden + log1p covariates | shared upstream | Methods §12.3 |
| **16** | **`16_aggregate_primary_sensitivity_v1/`** | **★ PRIMARY: aggregate Diff_any 4-test family + ASD-vs-SZ aggregate contrast + 4-layer robustness (standard/strict LOCO + top-N exclusion)** | **PRIMARY** | **Fig. 3a Block 1 + 3; §22, §26, §27** |
| 06 | `06_wgs_primary_L2/` | B′ logistic per-class (10 classes) + multinomial per-class ASD-vs-SZ heterogeneity | secondary localisation | Fig. 3b; §22 supporting; §27 secondary |
| 07 | `07_wgs_matched_static/` | Matched-static resampling specificity test (per-class) | sensitivity | Fig. 3c; §24 |
| 08 | `08_arraycgh_sample_burden/` | Japanese arrayCGH cohort burden | external | Fig. 3a Block 2; Fig. 3b external |
| 09 | `09_mssng_sample_burden/` | MSSNG ASD WGS GEE family-clustered burden | external | Fig. 3a Block 2; Fig. 3b external |
| 10 | `10_replication_2way_meta/` | Two-cohort fixed-effect IVW meta-analysis (arrayCGH + MSSNG, excluding WGS Discovery) | external | Fig. 3a Block 2 |
| 11 | `11_diff_all_vs_static/` | Diff_any vs Static per-bin specificity | sensitivity | Fig. 3d; §24 |
| 12 | `12_exon_exclusion_wgs/` | A1 + A2 exon-exclusion analyses in discovery WGS | sensitivity | §25 |
| 13 | `13_exon_exclusion_mssng/` | A1 + A2 exon-exclusion in MSSNG (Holm-adjusted) | sensitivity | §25 |
| 14 | `14_constraint_enrichment_v1/` | Three-cohort IVW constraint enrichment (Diff_any-bounded TAD-scale regulatory domains) | sensitivity (domain-level) | Fig. 3f; §29 |
| 99 | `99_verify_vs_draft/` | Numerical-reproducibility verification scripts | reproducibility | Methods §12 (transparency) |
| — | `common/` | Shared helper scripts (sample-info loading, covariate computation, BH-FDR utilities) | shared utility | — |

## Software dependencies

- Python 3.9+ (pandas, numpy, scipy.stats, statsmodels, scikit-learn for occasional helpers)
- R 4.2+ (readr, dplyr, tidyr, ggplot2, scales, gtable; for joint-layer also R-base BFGS multinomial implementation)
- BEDTools v2.30+
- SLURM workload manager (sbatch scripts)

## Input data and access

Reproducing the analyses additionally requires access to the cohort-specific input data and reference resources used in this study, including: WGS variant calls and structural-variant calls produced by the harmonised processing pipeline; ancestry principal components and depth metrics; the Japanese arrayCGH external cohort genotypes (Kushima et al.); MSSNG ASD WGS structural variants (de Rubeis lab; access via dbGaP/MSSNG portal); and the Heffel et al. snm3C-seq L2 differential boundary atlas (published bed files). All cohort-specific genotype data are subject to access controls described in the manuscript Data availability statement.

## Output

Each numbered step has its own `output_*/` subdirectory. The aggregate primary results (`16_aggregate_primary_sensitivity_v1/output/`) feed Supplementary Tables 5/6, Figure 3a–b, and Figure 5 (joint-layer ratios). The per-class secondary results (`06_wgs_primary_L2/output_v5/`) feed Supp Table 5 Panel A class-level rows.

## Endpoint mapping to figures and tables

- **Figure 3a Block 1** (4 primary tests): `16_aggregate_primary_sensitivity_v1/compute_aggregate_diff_any_4tests_v2.py`
- **Figure 3a Block 2** (External 2-cohort IVW): `08`, `09`, `10`
- **Figure 3a Block 3** (ASD-vs-SZ aggregate contrast): `16_aggregate_primary_sensitivity_v1/compute_aggregate_diff_any_asd_vs_sz_v2.py`
- **Figure 3b** (Secondary class-level localisation): `06_wgs_primary_L2/39_fit_B_prime_L2_and_specificity_v7.R` + external panel from `08`–`10`
- **Figure 3c** (Matched-static): `07_wgs_matched_static/tad_dynamic_boundary_specificity_test_v7.py`
- **Figure 3d** (Diff_any vs Static per-bin): `11_diff_all_vs_static/42_compute_diff_any_vs_static_v5.py`
- **Figure 3e** (Representative loci): exploratory locus illustration; selection script in `06_wgs_primary_L2/`
- **Figure 3f** (Constraint enrichment): `14_constraint_enrichment_v1/`
- **Supp Table 5 Panel A** (aggregate primary 4 tests + class-level): `16/` + `06/`
- **Supp Table 6 Panel A** (external IVW): `08`–`10`
- **Supp Table 6 Panel C** (Diff_any vs Static per-bin): `11`
- **Supp Table 6 Panel D** (exon-exclusion in WGS Discovery + MSSNG): `12`, `13` (arrayCGH excluded — see Sensitivity layers section above)
- **Supp Table 6 Panel E** (ASD-vs-SZ aggregate + class-level): `16/compute_aggregate_diff_any_asd_vs_sz_v2.py` + `06/40_fit_MNLogit_heterogeneity_v3.py` + arrayCGH cross-disorder (implemented within `08_arraycgh_sample_burden/tad_replication_arraycgh_v22.py` Phases 6B–6C)
- **Supp Table 6 Panel G** (constraint enrichment 3-cohort IVW): `14_constraint_enrichment_v1/`
- **Supp Table 6 Panel L** (4-layer aggregate signal robustness: standard/strict LOCO + top-N): `16/compute_loco_diff_any_sensitivity_v2.py` + `16/compute_strict_loco_top_exclusions_v1.py`
- **Figure 5** (joint-layer logistic regression): `joint_layer/asd_sz_layered_logistic_v6.R` (sibling repository directory)
- **Supp Table 11** (joint-layer): `joint_layer/asd_sz_layered_logistic_v6.R`

## Multiple-testing strategy

A hierarchical BH-FDR is used (Methods §12.4):
- **Aggregate primary** (4 tests; ASD/SZ × DEL/DUP): BH-FDR within 4-test family.
- **Class-level secondary** (10 classes per disorder per SV type): BH-FDR within 10-class family.
- **External replication** (9 discovery-positive classes): one-sided BH-FDR within 9-class family (the discovery-non-significant class and negative controls use two-sided tests).
- **Exon-exclusion** (Holm-adjusted one-sided P for MSSNG replication): family-wise error rate at the replication stage.

## Reproducibility verification

Run `99_verify_vs_draft/` after the pipeline completes to confirm that the primary L2 burden numbers and the aggregate Diff_any 4-test family numbers match the values reported in the manuscript and Supplementary Tables. Numerical mismatches indicate input-data or environment drift and should be reconciled before submission.

## Citation

If you use this code, please cite:

> Kushima I., et al. *Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia*. (in submission)
