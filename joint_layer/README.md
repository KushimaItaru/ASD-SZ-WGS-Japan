# Joint-layer integrative logistic regression analysis

This directory contains the integrative joint-layer logistic-regression analysis described in *Kushima et al., Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia*. The analysis evaluates whether the primary WGS-enabled signals (TAD Diff_any-boundary deletion burden and rare TRE outlier carriage) retain independent contributions when conditioned on established coding, CNV, polygenic, or repeat-expansion risk layers.

## Output mapping to manuscript

This directory generates the data underlying:
- **Figure 5a–b** — ASD and SZ joint-layer M0→Mfull progression OR forest plots
- **Figure 5c** — Full-model layer ORs across six progressive risk layers
- **Figure 5d** — Multinomial logistic regression with shared controls; ASD-vs-SZ TAD coefficient contrast
- **Figure 5e** — Static-boundary adjustment sensitivity (Diff_any-DEL alongside Static-DEL)
- **Supplementary Table 11** (Panels A–G) — full numerical results

## Script

**`asd_sz_layered_logistic_v6.R`** (single R script, ~31 KB)

Implementation summary:
- ASD + SZ unified joint-layer logistic regression (parallel binary models for each disorder)
- Six progressive risk layers added in nested order: M0 (B′-aligned base) → +TAD (Diff_any_DEL count) → +TRE (rare_any binary; cross-fitted EHdn or STRling) → +PRS (disorder-specific PRS-CSx, MHC-excluded, z-scored within Controls of joint subset) → +PTV (binary; rare PTV in disorder-specific risk-gene set) → +CNV (binary; rare exonic deletion in disorder-specific risk-gene set) → +NAHR-CNV (binary; carrier of disorder-specific NAHR-mediated genomic-disorder loci)
- Added-last likelihood-ratio test (LRT) for each layer's conditional contribution; per-layer P value reflects added-last (order-invariant) significance
- Multinomial logistic regression with shared controls (Controls = baseline; ASD = 1; SZ = 2) implemented in base R via maximum-likelihood estimation (BFGS optimizer with analytical gradient and Hessian inverse) for the formal Wald contrast on β_TAD,ASD − β_TAD,SZ; the Wald statistic is computed as Δβ / sqrt[var(β_ASD) + var(β_SZ) − 2·cov(β_ASD, β_SZ)] using the diagonal and off-diagonal vcov entries
- Static-boundary sensitivity model (log1p_Static_DEL added alongside Diff_any_DEL with the other five layers)
- glm convergence diagnostics; "possible_separation" warnings refer to ancillary covariates with extreme coefficients and do not affect inference on TAD or TRE main effects

The base model B′ covariates are: Sex_numeric + PC1–PC10 + log1p(total_del_bases) + log1p(total_gene_DEL); sequencing depth is additionally included for TRE-layer compatibility.

## Sample sizes

- **Joint subset for ASD-vs-Controls and SZ-vs-Controls binary models**: ASD 505, SZ 627, Controls 8,372 → N_total = 8,877 (ASD model) and 8,999 (SZ model)
- **Multinomial subset (shared controls; Fig. 5d)**: 8,372 Controls + 505 ASD + 627 SZ → N_total = 9,504
- The discrepancy in N between binary (8,877 / 8,999) and multinomial (9,504) reflects that the multinomial model uses all three classes simultaneously, whereas the binary models use only one disorder class plus controls

## Bug fix history (v5 → v6)

v5 used readr::read_tsv with default guess_max = 1000, which caused the Matched_SZ column class to be misinferred as logical when the first 1,000 rows were all NA (the first non-empty value occurred at row 1,153). This led to silent NA conversion of subsequent values and an erroneous count of 0 SZ × CNV carriers (true count: 39). v6 corrects this by setting guess_max = Inf for all relevant TSV reads, restoring correct typing and downstream layer counts.

## Output files (kept on NIG; not deposited in this repository)

The script produces nine TSV files plus a summary text file:
- `asd_sz_layered_v6.tsv` (Panel A; M0–Mfull progression, 28 rows)
- `asd_sz_layered_v6.added_last_lrt.tsv` (Panel B; 24 rows)
- `asd_sz_layered_v6.multinomial.tsv` (Panel C; 5 rows)
- `asd_sz_layered_v6.static_sensitivity.tsv` (Panel D; 4 rows)
- `asd_sz_layered_v6.carrier_counts.tsv` (Panel E)
- `asd_sz_layered_v6.cor_matrix.tsv` (Panel F; 84 rows; Pearson r among 7 layer variables in Controls)
- `asd_sz_layered_v6.diagnostics.tsv` (Panel G; 32 rows; glm convergence diagnostics)
- `asd_sz_layered_v6.summary.txt`
- `asd_sz_layered_v6.cohort.{ASD,SZ}.{EHdn,STRling}.tsv` (4 cohort-level cohort tables)

These output TSV files feed the corresponding panels of Supplementary Table 11 (xlsx; deposited as a Supplementary File with the manuscript). Output TSVs are not redistributed in this repository because they can be regenerated from the deposited script and the input data layers; the xlsx version with Cover and panel-level descriptions is the canonical published artifact.

## Key numerical results (verified against manuscript)

| Manuscript reference | Joint-layer output |
|---|---|
| ASD Mfull Diff_any-DEL OR (paper §33) | EHdn 1.392 (1.172–1.653); STRling 1.401 (1.180–1.665) |
| SZ Mfull Diff_any-DEL OR (paper §33) | EHdn 1.121 (0.938–1.340); STRling 1.118 (0.936–1.336) |
| ASD added-last LRT P, TAD layer (paper §33) | EHdn 4.20 × 10⁻⁵; STRling 2.99 × 10⁻⁵ |
| SZ added-last LRT P, TAD layer (paper §33) | EHdn 0.222; STRling 0.231 |
| Multinomial Wald P, β_TAD,ASD − β_TAD,SZ (paper §33) | 2.83 × 10⁻² (TAD-only caller-independent variant) |
| Δβ_TAD,ASD−SZ (paper §33) | 0.240 |
| Static-adjusted ASD Diff_any OR (paper §33) | EHdn 1.306; STRling 1.316 |
| Static-adjusted ASD Static OR (paper §33) | 1.010 (both callers) |
| Joint-subset N (paper §33) | 8,877 (ASD); 8,999 (SZ); 9,504 (multinomial) |

## Software dependencies

- R 4.2+ with packages: readr, dplyr, tidyr (no `nnet` package required; multinomial logistic regression is implemented in base R via `optim()` with the BFGS optimizer for portability across R installations)

## Input data and access

The script reads cohort-level layer files produced by upstream pipelines:
- TAD per-individual burden (from the TAD pipeline; tad/05_wgs_sample_burden output)
- TRE rare-outlier per-individual carriage (from cross-fitted EHdn and STRling pipelines)
- PRS scores (PRS-CSx; from the PRS pipeline)
- Rare PTV / exonic CNV / NAHR-CNV per-individual carrier indicators (from the rare-coding/CNV pipeline)

Cohort-specific genotype-derived inputs are subject to the access controls described in the manuscript Data availability statement.

## Citation

If you use this code, please cite:

> Kushima I., et al. *Shared rare tandem-repeat expansion burden and ASD-preferential developmental 3D-boundary deletion burden across autism and schizophrenia*. Manuscript in preparation, 2026.
