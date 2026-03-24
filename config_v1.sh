#!/bin/bash
# config_v1.sh
# - プロジェクト全体のパス・共通パラメータを集中管理
#
# 環境変数で上書き可能なパスにはデフォルト値を設定しています。
# 実行環境に合わせて export してから source してください。

set -euo pipefail

# ===== Project root =====
# このファイルは <PROJECT_ROOT>/config_v1.sh を想定
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ===== Common directories (created by setup) =====
export RES_DIR="${PROJECT_ROOT}/resources"
export WORK_DIR="${PROJECT_ROOT}/work"
export LOG_DIR="${PROJECT_ROOT}/logs"

export SAMPLE_LIST_DIR="${PROJECT_ROOT}/sample_lists"
export EHDN_OUT_DIR="${PROJECT_ROOT}/ehdn_output"
export MERGED_NOVEL_DIR="${PROJECT_ROOT}/merged_results_novel"
export ANALYSIS_NOVEL_DIR="${PROJECT_ROOT}/analysis_results_novel"

# ===== Input metadata =====
# SAMPLE_INFO: サンプルメタデータファイルへのパス（環境変数で上書き可能）
export SAMPLE_INFO="${SAMPLE_INFO:-${PROJECT_ROOT}/sampleInfo/GRIFIN_srWGS_SampleInfo.txt}"

# ===== CRAM base dirs (優先順に探索) =====
export CRAM_BASE_DIR1="${CRAM_BASE_DIR1:-/path/to/cram/case}"
export CRAM_BASE_DIR2="${CRAM_BASE_DIR2:-/path/to/cram/control}"

# ===== Reference =====
export REFERENCE_FASTA="${REFERENCE_FASTA:-/path/to/Homo_sapiens_assembly38.fasta}"

# ===== EHdn =====
export EHDN_BIN="${EHDN_BIN:-ExpansionHunterDenovo}"
export EHDN_MIN_ANCHOR_MAPQ="50"
export EHDN_MAX_IRR_MAPQ="40"
export EHDN_MIN_UNIT_LEN="3"
export EHDN_MAX_UNIT_LEN="8"

# ===== STRling（本burden解析では depth 計算のみ使用） =====
export DEPTH_DIR="${PROJECT_ROOT}/depth"
export DEPTH_PARTS_DIR="${DEPTH_DIR}/depth_parts"
export DEPTHS_ALL_TSV="${DEPTH_DIR}/depths_all.tsv"

# 深度の概算式の read length（srWGSで 151bp を想定）
export READ_LENGTH_BP="151"
export GENOME_SIZE_BP="3100000000"

# ===== Gene regions BED =====
# setupで resources/gene_regions_1kb_pad.bed を用意する想定
export GENE_REGIONS_BED="${RES_DIR}/gene_regions_1kb_pad.bed"

# ===== Normalization =====
export TARGET_DEPTH="40.0"

# ===== Slurm =====
export SLURM_PARTITION="${SLURM_PARTITION:-ncbn-cpu}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-ncbn-cpu}"

# EHdn job resources
export EHDN_CPUS="8"
export EHDN_MEM="32G"
export EHDN_TIME="48:00:00"
export EHDN_ARRAY_MAXPAR="40"

# Depth job resources
export DEPTH_CPUS="1"
export DEPTH_MEM="1G"
export DEPTH_TIME="1:00:00"
export DEPTH_ARRAY_MAXPAR="100"
