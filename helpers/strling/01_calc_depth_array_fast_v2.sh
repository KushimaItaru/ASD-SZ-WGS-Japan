#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
#SBATCH -J calc_depth_str1228
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-ARRAY_SIZE%100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 1:00:00
#SBATCH --mem=1G
#SBATCH --output=logs/depth/depth_%A_%a.out
#SBATCH --error=logs/depth/depth_%A_%a.err

# strling/01_calc_depth_array_fast_v2.sh (FIXED)
# - sample_lists/ehdn_all_samples.tsv を対象に深度計算（family_member含む）
# - Slurm実行環境でも config_v1.sh を確実に読むため、PROJECT_ROOT を絶対パスで指定
# - 出力: depth/depth_parts/<SampleID>.tsv
# - 実行時間を記録

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

# ★重要：Slurm上で確実に見つかるようプロジェクトルートを絶対パスで固定
PROJECT_ROOT="/lustre12/home/kushima-pg/str_12282025"  # CONFIGURE
source "${PROJECT_ROOT}/config_v1.sh"

mkdir -p "${LOG_DIR}/depth" "${DEPTH_PARTS_DIR}"

SAMPLE_LIST="${SAMPLE_LIST_DIR}/ehdn_all_samples.tsv"

LINE=$((SLURM_ARRAY_TASK_ID + 1))
REC=$(sed -n "${LINE}p" "${SAMPLE_LIST}" || true)
if [ -z "${REC}" ]; then
  exit 0
fi

SAMPLE_ID=$(echo "${REC}" | cut -f1)
CRAM_PATH=$(echo "${REC}" | cut -f3)
OUT="${DEPTH_PARTS_DIR}/${SAMPLE_ID}.tsv"

# 既にあるならスキップ
if [ -s "${OUT}" ]; then
  exit 0
fi

if [ ! -f "${CRAM_PATH}" ]; then
  echo -e "${SAMPLE_ID}\t0.0" > "${OUT}"
  exit 0
fi

samtools idxstats "${CRAM_PATH}" | \
awk -v id="${SAMPLE_ID}" -v rl="${READ_LENGTH_BP}" -v gs="${GENOME_SIZE_BP}" '
  {sum += $3}
  END{
    depth = (sum * rl) / gs;
    printf "%s\t%.4f\n", id, depth
  }
' > "${OUT}"

END=$(date +%s)
echo "[$(TS)] [DONE] ${SAMPLE_ID} depth computed. Elapsed=$((END-START))s"
