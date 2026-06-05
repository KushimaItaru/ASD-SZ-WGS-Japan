# Code Availability — Draft text for manuscript

This file contains the canonical Code availability paragraph used in the manuscript (matching paper v313), together with shorter alternatives for length-constrained venues.

## Recommended version (matches manuscript Code availability paragraph)

Code supporting all primary and supporting analyses is available at
https://github.com/KushimaItaru/ASD-SZ-WGS-Japan. The repository contains,
at the time of submission, the modules covering: (i) the TRE cross-fitting
and rare-burden pipeline (`tre/`; entry-point wrapper scripts for
ExpansionHunter Denovo and STRling, together with cross-fitting,
rare-frequency thresholding, and burden-regression helper scripts);
(ii) the TAD boundary disruption analysis with a two-tier endpoint
structure: an architecture-level primary endpoint
(`tad/16_aggregate_primary_sensitivity_v1`; aggregate Diff_any boundary-bin
four-test family ASD/SZ × DEL/DUP, the aggregate ASD-versus-SZ contrast,
and four-layer aggregate-signal robustness comprising standard and strict
leave-one-class-out plus top-N high-burden sample and event exclusion), a
shared upstream pipeline constructing the boundary universe and
per-individual burden (`tad/01_heffel_boundary_master` through
`tad/05_wgs_sample_burden`), and sensitivity layers comprising
matched-static resampling with chromosome × gene-density × segmental-duplication
matching (`tad/07_wgs_matched_static`), Diff_any-versus-Static per-bin
specificity testing (`tad/11_diff_all_vs_static`), and exon-exclusion
two-test analyses across discovery WGS and MSSNG
(`tad/12_exon_exclusion_wgs`, `tad/13_exon_exclusion_mssng`);
(iii) the external replication and IVW meta-analysis pipeline for the
arrayCGH and MSSNG cohorts (`tad/08_arraycgh_sample_burden`,
`tad/09_mssng_sample_burden`, `tad/10_replication_2way_meta`; within-array
IVW, GEE family-clustered logistic regression, fixed-effect three-cohort
IVW combining WGS Discovery, arrayCGH, and MSSNG, and the Diff_any-bounded
TAD-scale regulatory-domain constraint-enrichment analysis in
`tad/14_constraint_enrichment_v1`); (iv) the secondary class-level
localisation (`tad/06_wgs_primary_L2`; per-class B′ logistic regression and
per-class multinomial heterogeneity testing with shared-control covariance
and correlation-adjusted Stouffer–Brown class-level global directional
testing); and (v) numerical-reproducibility verification scripts that
re-derive the primary L2 burden numbers reported in the manuscript
(`tad/99_verify_vs_draft`). The integrative joint-layer logistic regression
for ASD and SZ is available at
`joint_layer/asd_sz_layered_logistic_v6.R`. Reproducing the full analyses
additionally requires access to the cohort-specific input data and
reference resources used in this study.

## Shorter version (if space is limited)

Code for the tandem repeat expansion (TRE) pipeline, the TAD boundary
disruption analysis (aggregate Diff_any primary endpoint, shared upstream
pipeline, sensitivity layers, external replication and three-cohort IVW
meta-analysis, secondary class-level localisation, and constraint
enrichment), and the integrative joint-layer logistic regression for ASD
and SZ is available at https://github.com/KushimaItaru/ASD-SZ-WGS-Japan.
The repository contains entry-point wrapper scripts, helper scripts,
configuration files, and module-level documentation needed to reproduce
the main analyses. Reproducing the full analyses additionally requires
access to the cohort-specific input data and reference resources used in
this study.

## Minimal version (very space-constrained)

Code is available at https://github.com/KushimaItaru/ASD-SZ-WGS-Japan,
covering the TRE rare-burden pipeline (`tre/`), the TAD boundary disruption
pipeline with primary aggregate-endpoint, sensitivity layers, external
replication, and class-level localisation modules (`tad/`), and the
joint-layer integrative logistic regression
(`joint_layer/asd_sz_layered_logistic_v6.R`). Reproducing the full
analyses additionally requires access to the cohort-specific input data
and reference resources used in this study.

---

## GitHub About field (repository description)

Publication-ready analysis pipelines for ASD/SZ WGS: rare tandem-repeat
expansion (TRE) burden, cell-type-resolved TAD boundary disruption with
two-tier endpoint structure, and joint-layer integrative logistic
regression.

## GitHub README opening paragraph

This repository contains the publication-ready workflows underlying
Kushima et al. (2026 (manuscript under review)): (i) the rare tandem-repeat
expansion (TRE) burden pipeline based on ExpansionHunter Denovo and
STRling with 5-fold cross-fitted regression; (ii) the cell-type-resolved
TAD boundary disruption pipeline with an architecture-level aggregate
primary endpoint, secondary class-level localisation, sensitivity layers,
external replication and three-cohort IVW meta-analysis, and constraint
enrichment; and (iii) the integrative joint-layer logistic regression
spanning six progressive risk layers for ASD and SZ. Each module ships
with entry-point wrappers, helper scripts, configuration, and
documentation sufficient to reproduce the main analyses given access to
the cohort-specific input data and reference resources.
