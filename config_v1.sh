#!/bin/bash
# config_v1.sh
# - プロジェクト全体のパス・共通パラメータを集中管理

set -euo pipefail

# ===== Project root =====
# このファイルは ~/str_12282025/config_v1.sh を想定
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
export SAMPLE_INFO="/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt"

# ===== CRAM base dirs (優先順に探索) =====
export CRAM_BASE_DIR1="/lustre12/home/grifinpd-pg/analysis/parabricks"
export CRAM_BASE_DIR2="/lustre12/home/ncbn-share-pg/control_genome/pb3.1.0/results"

# ===== Reference =====
export REFERENCE_FASTA="/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta"

# ===== EHdn =====
export EHDN_BIN="/home/kushima-pg/ehdn_analysis/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo"
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

# ===== Slurm (NCBN) =====
export SLURM_PARTITION="ncbn-cpu"
export SLURM_ACCOUNT="ncbn-cpu"

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
