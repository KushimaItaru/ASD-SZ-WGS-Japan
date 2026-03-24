#!/usr/bin/env bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# strling/00_strling_genomewide_config_v1.sh
# - Centralized management of output directories, references, and sample lists for STRling (genome-wide)
# - Default sample list set to ehdn_all_samples.tsv (including family_member)
# - Define runtime environment (conda/strling) and chromosome list for merge

set -euo pipefail

# ---- Project root (fixed)----
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "${PROJECT_ROOT}/config_v1.sh"

# ---- STRling env / binary ----
export STRLING_USE_CONDA=true
export STRLING_CONDA_ENV="strling_env"
export STRLING_BIN="${STRLING_BIN:-strling}"

# ---- Input sample list ----
# ★Default: include family_member (same as EHdn)
export STRLING_SAMPLES_TSV="${SAMPLE_LIST_DIR}/ehdn_all_samples.tsv"

# Use below for case-control only (manual switch)
# export STRLING_SAMPLES_TSV="${SAMPLE_LIST_DIR}/casecontrol_samples.tsv"

# ---- Output dirs (within this project) ----
export STRLING_OUT_ROOT="${STRLING_OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
export STRLING_BINS_DIR="${STRLING_BINS_DIR:-${STRLING_OUT_ROOT}/bins}"
export STRLING_RES_DIR="${STRLING_RES_DIR:-${STRLING_OUT_ROOT}/str-results}"
export STRLING_LOG_DIR="${STRLING_LOG_DIR:-${STRLING_OUT_ROOT}/logs}"

mkdir -p "${STRLING_BINS_DIR}" "${STRLING_RES_DIR}" "${STRLING_LOG_DIR}"

# ---- Reference FASTA ----
export STRLING_REFERENCE_FASTA="${REFERENCE_FASTA}"

# ---- STRling parameters (reasonable defaults) ----
export STRLING_PROPORTION_REPEAT="0.8"
export STRLING_MIN_MAPQ="40"

# merge parameters
export STRLING_MIN_SUPPORT="5"
export STRLING_MIN_CLIP_TOTAL="0"
export STRLING_WINDOW="-1"

# Chromosomes (hg38 chr-prefix)
export STRLING_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
                       chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
                       chr21 chr22 chrX chrY)

# ---- Slurm resources ----
# extract
export STRLING_EXTRACT_CPUS="1"
export STRLING_EXTRACT_MEM="32G"
export STRLING_EXTRACT_TIME="12:00:00"
export STRLING_EXTRACT_MAXPAR="40"

# merge
export STRLING_MERGE_CPUS="1"
export STRLING_MERGE_MEM="128G"
export STRLING_MERGE_TIME="24:00:00"
export STRLING_MERGE_MAXPAR="24"   # per chr, array=1-24

# ---- Downstream alias (config → downstream variable names) ----
# Downstream scripts reference variable names without STRLING_ prefix,
# so export aliases for compatibility.
export OUT_ROOT="${STRLING_OUT_ROOT}"
export BINS_DIR="${STRLING_BINS_DIR}"
export STR_RES_DIR="${STRLING_RES_DIR}"
export LOG_DIR="${STRLING_LOG_DIR}"
export CHROMOSOMES=("${STRLING_CHROMS[@]}")
export WINDOW="${STRLING_WINDOW}"
export MIN_SUPPORT="${STRLING_MIN_SUPPORT}"
export MIN_CLIP_TOTAL="${STRLING_MIN_CLIP_TOTAL}"
export MIN_MAPQ="${STRLING_MIN_MAPQ}"
export PROPORTION_REPEAT="${STRLING_PROPORTION_REPEAT}"
