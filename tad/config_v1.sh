#!/bin/bash
# config_v1.sh
# - TAD pipeline (tad04292026) の中央 path/parameter 設定
# - 各 path は env var 優先で、設定なしの場合は NIG default を使用
# - TRE module の config_v1.sh と同じ ${VAR:-default} pattern
#
# Usage:
#   - 何もしなくても NIG 環境では default が effective
#   - 外部 reproducer は事前に export VAR=... してから sbatch するか、
#     `source config_v1.sh` で env var を current shell に exposed する
#
# 関連: common/paths_v1.py が同じ env var を読み出す

set -euo pipefail

# ============================================================================
# Project root (auto-detect; this file at <PROJECT_ROOT>/config_v1.sh)
# ============================================================================
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# Pipeline data root (where module output dirs live; defaults to NIG path)
# ============================================================================
# Python paths_v1.py が PIPELINE_ROOT env var を読む
export PIPELINE_ROOT="${PIPELINE_ROOT:-/lustre12/home/kushima-pg/tad04292026}"  # CONFIGURE

# ============================================================================
# External raw data resources (Heffel atlas, sample metadata, references)
# ============================================================================
export HEFFEL_ROOT="${HEFFEL_ROOT:-/home/kushima-pg/resource/heffel2024nature}"  # CONFIGURE
export SAMPLE_INFO="${SAMPLE_INFO:-/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt}"  # CONFIGURE
export REFERENCE_GENOME="${REFERENCE_GENOME:-/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta}"  # CONFIGURE
export GENCODE_GTF="${GENCODE_GTF:-/lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz}"  # CONFIGURE

# ============================================================================
# CRAM file path templates (for accessing per-sample CRAM)
# ============================================================================
export CRAM_TEMPLATE_GRIFIN="${CRAM_TEMPLATE_GRIFIN:-/lustre12/home/grifinpd-pg/analysis/parabricks/{sample_id}/{sample_id}.cram}"  # CONFIGURE
export CRAM_TEMPLATE_NCBN="${CRAM_TEMPLATE_NCBN:-/lustre12/home/ncbn-share-pg/control_genome/pb3.1.0/results/{sample_id}/{sample_id}.cram}"  # CONFIGURE

# ============================================================================
# External replication cohorts (arrayCGH, MSSNG)
# ============================================================================
export ARRAYCGH_DATA_ROOT="${ARRAYCGH_DATA_ROOT:-/lustre12/home/kushima-pg/arraycgh_data}"  # CONFIGURE
export MSSNG_DATA_ROOT="${MSSNG_DATA_ROOT:-/lustre12/home/kushima-pg/mssng_data}"  # CONFIGURE

# ============================================================================
# CNV caller pipeline outputs (AnnotSV-derived)
# ============================================================================
export CNV_ANALYSIS_ROOT="${CNV_ANALYSIS_ROOT:-/lustre12/home/kushima-pg/cnv_01012026/analysis}"  # CONFIGURE
export CURATED_GD_FILE="${CURATED_GD_FILE:-/lustre12/home/kushima-pg/cnv_01012026/curated_genomic_disorder_cnv_loci_v3.txt}"  # CONFIGURE

# ============================================================================
# Genome annotation resources
# ============================================================================
export SEGDUP_BED="${SEGDUP_BED:-/lustre12/home/kushima-pg/annotationInfo/segdup_hg38_sorted_merged.bed}"  # CONFIGURE
export EXCLUSION_BED="${EXCLUSION_BED:-/lustre12/home/kushima-pg/resource/hg38_cnv_exclusion_regions.bed}"  # CONFIGURE
export PCA_EIGENVEC="${PCA_EIGENVEC:-/home/kushima-pg/PRS/population_stratfication_09012025/results_popstrat_20250903_v7/pca_jpn/pca.eigenvec}"  # CONFIGURE
export LINEAGE_XLSX="${LINEAGE_XLSX:-/lustre12/home/kushima-pg/heffel_deep_analysis_03242026/lineage_stage_clusters_v3.xlsx}"  # CONFIGURE

# ============================================================================
# Pipeline parameters
# ============================================================================
# Bin size (Heffel atlas raw resolution, 25 kb)
export BIN_SIZE_BP="${BIN_SIZE_BP:-25000}"

# Top-1% sample-level CNV count QC threshold (Module 04 v10 / Module 03 v8+)
export CNV_TOP_PCT_QUANTILE="${CNV_TOP_PCT_QUANTILE:-0.99}"

# B' logistic regression min cell count (sparse exposure guard)
export MIN_CELL_COUNT="${MIN_CELL_COUNT:-5}"

# matched-static resampling iterations (Module 07)
export MATCHED_STATIC_RESAMPLES="${MATCHED_STATIC_RESAMPLES:-10000}"

# ============================================================================
# SLURM defaults (NIG)
# ============================================================================
export SLURM_PARTITION="${SLURM_PARTITION:-ncbn-cpu}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-ncbn-cpu}"

# ============================================================================
# Conda environment name (for analytical scripts)
# ============================================================================
export CONDA_ENV_NAME="${CONDA_ENV_NAME:-tad-burden}"

# ============================================================================
# (Optional) Print active configuration when sourced verbosely
# ============================================================================
if [[ "${TAD_CONFIG_VERBOSE:-0}" == "1" ]]; then
    echo "[config_v1.sh] PROJECT_ROOT=${PROJECT_ROOT}"
    echo "[config_v1.sh] PIPELINE_ROOT=${PIPELINE_ROOT}"
    echo "[config_v1.sh] HEFFEL_ROOT=${HEFFEL_ROOT}"
    echo "[config_v1.sh] SAMPLE_INFO=${SAMPLE_INFO}"
    echo "[config_v1.sh] REFERENCE_GENOME=${REFERENCE_GENOME}"
fi
