"""
tad04292026/common/paths_v1.py

Centralized path constants for the TAD boundary re-analysis pipeline (paper-aligned
public version).

Usage in each script:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from common.paths_v1 import *
    # then inside main(): ensure_output_dirs()

History:
  2026-04-21: BIN_SIZE_BP 40_000 -> 25_000 (Heffel raw data is 25kb bins).
              Added export_r_manifest() for R-side single-source-of-truth.
              ensure_output_dirs() is now intended to be called from each
              script's main(), not at import time.
  2026-04-29: Path abstraction (TRE pattern: os.environ.get(VAR, NIG_default)).
              All hardcoded paths wrapped in env-var override; default values
              remain the actual NIG paths so the pipeline runs unchanged on the
              supercomputer with no env vars set, while external reproducers can
              override via env vars (e.g. `export PIPELINE_ROOT=/their/path`).
              Each external path is marked with `# CONFIGURE` for clarity.
              Module 02 v3 output (S2-S5 sensitivity grouping removed) is now
              the default: OUT_02_BIN_L2_ANNOT / bin_l2_annotation_v3.tsv.gz.
              PIPELINE_ROOT default updated tad04292026 -> tad04292026.
"""

from __future__ import annotations
import json
import os
from pathlib import Path


# ============================================================================
# Helper: env var with NIG default
# ============================================================================
def _envpath(var: str, default: str) -> Path:
    """Return Path(os.environ.get(var, default)). All external paths use this."""
    return Path(os.environ.get(var, default))


# ============================================================================
# RAW DATA ROOTS (external, read-only)
# ============================================================================

HEFFEL_ROOT = _envpath(
    "HEFFEL_ROOT",
    "/home/kushima-pg/resource/heffel2024nature",
)  # CONFIGURE
HEFFEL_DOMAIN_BOUNDARIES = HEFFEL_ROOT / "domain_boundaries"
HEFFEL_L2_DIFF_DOMAIN_BOUNDARIES = HEFFEL_ROOT / "L2_diff_domain_boundaries"

SAMPLE_INFO = _envpath(
    "SAMPLE_INFO",
    "/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt",
)  # CONFIGURE

REFERENCE_GENOME = _envpath(
    "REFERENCE_GENOME",
    "/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta",
)  # CONFIGURE
REFERENCE_GENOME_FAI = Path(str(REFERENCE_GENOME) + ".fai")

CRAM_TEMPLATE_GRIFIN = os.environ.get(
    "CRAM_TEMPLATE_GRIFIN",
    "/lustre12/home/grifinpd-pg/analysis/parabricks/{sample_id}/{sample_id}.cram",
)  # CONFIGURE
CRAM_TEMPLATE_NCBN = os.environ.get(
    "CRAM_TEMPLATE_NCBN",
    "/lustre12/home/ncbn-share-pg/control_genome/pb3.1.0/results/{sample_id}/{sample_id}.cram",
)  # CONFIGURE

ARRAYCGH_DATA_ROOT = _envpath(
    "ARRAYCGH_DATA_ROOT",
    "/lustre12/home/kushima-pg/arraycgh_data",
)  # CONFIGURE
MSSNG_DATA_ROOT = _envpath(
    "MSSNG_DATA_ROOT",
    "/lustre12/home/kushima-pg/mssng_data",
)  # CONFIGURE

GENCODE_GTF = _envpath(
    "GENCODE_GTF",
    "/lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz",
)  # CONFIGURE


# --- CNV caller pipeline outputs (AnnotSV etc.) ---
CNV_ANALYSIS_ROOT = _envpath(
    "CNV_ANALYSIS_ROOT",
    "/lustre12/home/kushima-pg/cnv_01012026/analysis",
)  # CONFIGURE
ANNOTSV_TABLE     = CNV_ANALYSIS_ROOT / "annotsv_results" / "annotsv_annotated_autosomal.tsv"
COMMON_CNV_TABLE  = CNV_ANALYSIS_ROOT / "annotsv_results" / "common_cnv_sites_0.1pct.tsv"
CNV_SAMPLE_COUNTS = CNV_ANALYSIS_ROOT / "cnv_sample_counts.tsv"
CURATED_GD_FILE = _envpath(
    "CURATED_GD_FILE",
    "/lustre12/home/kushima-pg/cnv_01012026/curated_genomic_disorder_cnv_loci_v3.txt",
)  # CONFIGURE

# --- Genome annotation resources ---
SEGDUP_BED = _envpath(
    "SEGDUP_BED",
    "/lustre12/home/kushima-pg/annotationInfo/segdup_hg38_sorted_merged.bed",
)  # CONFIGURE
EXCLUSION_BED = _envpath(
    "EXCLUSION_BED",
    "/lustre12/home/kushima-pg/resource/hg38_cnv_exclusion_regions.bed",
)  # CONFIGURE
PCA_EIGENVEC = _envpath(
    "PCA_EIGENVEC",
    "/home/kushima-pg/PRS/population_stratfication_09012025/results_popstrat_20250903_v7/pca_jpn/pca.eigenvec",
)  # CONFIGURE

LINEAGE_XLSX = _envpath(
    "LINEAGE_XLSX",
    "/lustre12/home/kushima-pg/heffel_deep_analysis_03242026/lineage_stage_clusters_v3.xlsx",
)  # CONFIGURE

# ============================================================================
# PIPELINE ROOT (this re-analysis tree)
# ============================================================================
# 2026-04-29: PIPELINE_ROOT default updated tad04292026 -> tad04292026
# (paper-aligned working repo with S2-S5 grouping removed and other paper-irrelevant
# processes pruned).

PIPELINE_ROOT = _envpath(
    "PIPELINE_ROOT",
    "/lustre12/home/kushima-pg/tad04292026",
)  # CONFIGURE

OUT_01_BOUNDARY_MASTER         = PIPELINE_ROOT / "01_heffel_boundary_master/output_v9"
# 2026-04-29: Module 02 v2 -> v3 (S2-S5 sensitivity grouping removed).
OUT_02_BIN_L2_ANNOT            = PIPELINE_ROOT / "02_bin_l2_annotation/output_v3"
OUT_03_WGS_SV_EVENTS           = PIPELINE_ROOT / "03_wgs_sv_events/output_v9"
OUT_04_SV_BOUNDARY_OVERLAP     = PIPELINE_ROOT / "04_wgs_sv_boundary_overlap/output_v10"
# 2026-04-29: Module 05 v3 -> v4, Module 06 v4 -> v5 (S2-S5 sensitivity grouping removed).
# Note: Pattern-specific output (Pattern B/C) is constructed in scripts from PATTERN env var:
#       Module 05 Pattern B/C: output_v6_{patB,patC}
#       Module 06 Pattern B/C: output_v7_{patB,patC}
OUT_05_SAMPLE_BURDEN           = PIPELINE_ROOT / "05_wgs_sample_burden/output_v4"  # primary (Pattern A, v4)
OUT_06_B_PRIME_L2              = PIPELINE_ROOT / "06_wgs_primary_L2/output_v5"     # primary (Pattern A, v5)
OUT_07_MATCHED_STATIC          = PIPELINE_ROOT / "07_wgs_matched_static/output_v6"
OUT_08_ARRAYCGH_BURDEN         = PIPELINE_ROOT / "08_arraycgh_sample_burden/output_v22"
OUT_09_MSSNG_BURDEN            = PIPELINE_ROOT / "09_mssng_sample_burden/output_v18"
OUT_09_MSSNG_META              = PIPELINE_ROOT / "09_mssng_sample_burden/meta_v12"
OUT_10_REPLICATION_META        = PIPELINE_ROOT / "10_replication_2way_meta/output_v5"
OUT_11_DIFF_ALL_VS_STATIC      = PIPELINE_ROOT / "11_diff_all_vs_static/output_v5"
OUT_12_WGS_EXONFREE            = PIPELINE_ROOT / "12_exon_exclusion_wgs/output_v3"
OUT_13_MSSNG_EXONFREE          = PIPELINE_ROOT / "13_exon_exclusion_mssng/output_v3"
OUT_14_CONSTRAINT_ENRICHMENT   = PIPELINE_ROOT / "14_constraint_enrichment_v1/output_v1"
OUT_15_L2_JACCARD              = PIPELINE_ROOT / "15_l2_jaccard_v1/output"
OUT_99_VERIFY                  = PIPELINE_ROOT / "99_verify_vs_draft/output_v1"

# Step 01
F_01_BOUNDARY_MASTER       = OUT_01_BOUNDARY_MASTER / "heffel_boundary_master_v9.tsv.gz"
F_01_BOUNDARY_MASTER_BED   = OUT_01_BOUNDARY_MASTER / "heffel_boundary_master_v9.bed"
F_01_BIN_PROPERTIES        = OUT_01_BOUNDARY_MASTER / "bin_properties_v9.tsv.gz"

# Step 02 (v3: S2-S5 removed; group_primary only)
F_02_BIN_L2_ANNOTATION     = OUT_02_BIN_L2_ANNOT / "bin_l2_annotation_v3.tsv.gz"

# Step 03 (Pattern A primary; Pattern B/C constructed in scripts)
F_03_WGS_SV_EVENTS         = OUT_03_WGS_SV_EVENTS / "wgs_rare_sv_events_v9_patA.tsv.gz"
F_03_WGS_SV_EVENTS_SUMMARY = OUT_03_WGS_SV_EVENTS / "wgs_rare_sv_events_v9_patA_summary.tsv"

# Step 04 (Pattern A primary)
F_04_EVENT_OVERLAP           = OUT_04_SV_BOUNDARY_OVERLAP / "sample_boundary_event_overlap_v10.tsv.gz"
F_04_SAMPLE_BURDEN_COV       = OUT_04_SV_BOUNDARY_OVERLAP / "sample_mechanistic_burden_v10.tsv"
F_04_SAMPLE_BURDEN_COV_SUM   = OUT_04_SV_BOUNDARY_OVERLAP / "sample_mechanistic_burden_summary_v10.tsv"

# Step 05 (Pattern A primary, v4 = S2-S5 削除済; Pattern B/C constructed in scripts)
F_05_SAMPLE_BURDEN_L2      = OUT_05_SAMPLE_BURDEN / "sample_burden_L2_and_specificity_v4.tsv"
F_05_SAMPLE_BURDEN_L2_SUM  = OUT_05_SAMPLE_BURDEN / "sample_burden_L2_and_specificity_summary_v4.tsv"

# Step 06 (Pattern A primary, v5 = S2-S5 削除済)
F_06_B_PRIME_L2_RESULTS    = OUT_06_B_PRIME_L2 / "B_prime_L2_classes_results_v5.tsv"
F_06_COVARIATES            = OUT_06_B_PRIME_L2 / "covariates_v5.tsv.gz"

# Step 07 (matched-static)
F_07_MATCHED_STATIC_MAIN   = OUT_07_MATCHED_STATIC / "matched_static_results_v6.tsv"

# Step 08 (arrayCGH)
F_08_ARRAYCGH_ASD          = OUT_08_ARRAYCGH_BURDEN / "tad_asd_vs_cont_v22.tsv"
F_08_ARRAYCGH_SZ           = OUT_08_ARRAYCGH_BURDEN / "tad_sz_vs_cont_v22.tsv"

# Step 09 (MSSNG)
F_09_MSSNG_BURDEN          = OUT_09_MSSNG_BURDEN / "mssng_sample_burden_v18.tsv.gz"
F_09_MSSNG_META            = OUT_09_MSSNG_META / "meta_results_ivw_v12.tsv"

# Step 10 (3-cohort IVW meta)
F_10_META_3WAY             = OUT_10_REPLICATION_META / "meta_results_ivw_v5.tsv"

# Step 11 (Diff_any vs Static; renamed from diff_all -> diff_any in v5)
F_11_DIFF_ANY_VS_STATIC    = OUT_11_DIFF_ALL_VS_STATIC / "meta_diff_any_vs_static_v5.tsv"

# Step 12 (WGS exon-free)
F_12_WGS_EXONFREE          = OUT_12_WGS_EXONFREE / "wgs_exon_exclusion_v3.tsv"

# Step 13 (MSSNG exon-free)
F_13_MSSNG_EXONFREE        = OUT_13_MSSNG_EXONFREE / "mssng_exon_exclusion_v3.tsv"

# Step 99 (verification)
F_99_PIPELINE_CONSISTENCY_REPORT = OUT_99_VERIFY / "pipeline_consistency_v1_report.tsv"
F_99_PIPELINE_CONSISTENCY_LOG    = OUT_99_VERIFY / "pipeline_consistency_v1.log"
F_99_PAPER_NUMBERS_REPORT  = OUT_99_VERIFY / "paper_numbers_verification_v1.tsv"
F_99_PAPER_NUMBERS_LOG     = OUT_99_VERIFY / "paper_numbers_verification_v1.log"

# R-side single-source-of-truth manifest
R_MANIFEST_JSON = PIPELINE_ROOT / "common" / "paths_v1_manifest.json"

# ============================================================================
# OUTPUT-DIR HELPER (call explicitly from each script's main(), not on import)
# ============================================================================

_OUTPUT_DIRS = [
    OUT_01_BOUNDARY_MASTER, OUT_02_BIN_L2_ANNOT, OUT_03_WGS_SV_EVENTS,
    OUT_04_SV_BOUNDARY_OVERLAP, OUT_05_SAMPLE_BURDEN, OUT_06_B_PRIME_L2,
    OUT_07_MATCHED_STATIC, OUT_08_ARRAYCGH_BURDEN, OUT_09_MSSNG_BURDEN,
    OUT_09_MSSNG_META, OUT_10_REPLICATION_META, OUT_11_DIFF_ALL_VS_STATIC,
    OUT_12_WGS_EXONFREE, OUT_13_MSSNG_EXONFREE, OUT_14_CONSTRAINT_ENRICHMENT,
    OUT_15_L2_JACCARD, OUT_99_VERIFY,
]

def ensure_output_dirs() -> None:
    for d in _OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# ============================================================================
# R MANIFEST (single source of truth for R scripts)
# ============================================================================

def export_r_manifest(manifest_path: Path | None = None) -> Path:
    """Write a JSON manifest containing all path constants consumed by R scripts.

    R side (e.g. 39_fit_B_prime_L2_and_specificity_v7.R) loads this via
      paths <- jsonlite::fromJSON("paths_v1_manifest.json")
      DEFAULT_BURDEN <- paths$F_05_SAMPLE_BURDEN_L2
    to avoid hand-maintained path duplication between Python and R.
    """
    if manifest_path is None:
        manifest_path = R_MANIFEST_JSON
    manifest = {
        # input files R needs
        "F_05_SAMPLE_BURDEN_L2": str(F_05_SAMPLE_BURDEN_L2),
        "F_05_SAMPLE_BURDEN_L2_SUM": str(F_05_SAMPLE_BURDEN_L2_SUM),
        # output locations R writes
        "OUT_06_B_PRIME_L2": str(OUT_06_B_PRIME_L2),
        "F_06_B_PRIME_L2_RESULTS": str(F_06_B_PRIME_L2_RESULTS),
        "F_06_COVARIATES": str(F_06_COVARIATES),
        # constants
        "BIN_SIZE_BP": BIN_SIZE_BP,
        "L2_CLASSES": L2_CLASSES,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return manifest_path

# ============================================================================
# CONSTANTS
# ============================================================================

# CamelCase-dash names (Heffel Nature 2024 L2 labels). 10 classes after excluding HPC_Astro.
L2_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]

B_PRIME_COVARIATES = [
    "Sex",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PC6", "PC7", "PC8", "PC9", "PC10",
    "log1p_total_del_bases",
    "log1p_total_gene_DEL",
]

# 2026-04-21: corrected from 40_000 to 25_000 to match Heffel raw data
#             (HPC_Exc-CA_diffbound.bed.gz bins are 25 kb; h5ad var coordinates confirmed 25 kb).
BIN_SIZE_BP = 25_000

# Make common/ a package
